module fdot_simplex
!! Network simplex, which assigns integer ranks to the nodes of a directed acyclic graph so that each
!! edge is at least `minlen` long, and the sum of `weight` times length over the edges is smallest.
!! This is a port of `lib/common/ns.c` of Graphviz, which is under the Eclipse Public License 2.0.
!!
!! The solver works on its own graph in plain arrays, so that the ranking and the positioning stages
!! can both use it. Searches that are recursive in the first versions of Graphviz use explicit stacks
!! here, so that a large graph cannot overflow the stack of a thread.
use iso_fortran_env, only: int64
implicit none
private

public :: ns_t, ns_allocate, ns_finish_lists, ns_solve

integer, parameter :: SEARCHSIZE = 30   !! how many negative edges to look at before choosing one

type :: ns_t
    !! A graph for the solver, and what the solver keeps while it works. Nodes and edges are numbered
    !! from 1, and 0 stands for none.
    integer :: n = 0                        !! number of nodes
    integer :: m = 0                        !! number of edges
    integer, allocatable :: tail(:)         !! where each edge starts
    integer, allocatable :: head(:)         !! where each edge ends
    integer, allocatable :: minlen(:)       !! the smallest length of each edge
    integer, allocatable :: weight(:)       !! the weight of each edge
    integer, allocatable :: rank(:)         !! the rank of each node: the input and the result
    logical, allocatable :: normal(:)       !! false for a node that is not a node of the drawing
    integer, allocatable :: out_start(:)    !! where the out-edges of each node start in `out_list`
    integer, allocatable :: in_start(:)     !! where the in-edges of each node start in `in_list`
    integer, allocatable :: out_list(:)     !! the out-edges of the nodes
    integer, allocatable :: in_list(:)      !! the in-edges of the nodes

    integer(int64), allocatable, private :: cutvalue(:)
    integer, allocatable, private :: tree_index(:), tree_edge(:), tree_out(:), tree_in(:)
    integer, allocatable, private :: tree_out_size(:), tree_in_size(:)
    integer, allocatable, private :: par(:), low(:), lim(:), priority(:), subtree(:)
    integer, allocatable, private :: st_rep(:), st_size(:), st_heap(:), st_par(:), heap(:)
    integer, allocatable, private :: stack(:, :)
    integer, private :: tree_count = 0, search_index = 0, search_size = SEARCHSIZE, heap_size = 0
end type ns_t

contains

!> Make room for a graph. The caller then fills `tail`, `head`, `minlen`, `weight` and `rank`, and
!> `out_start`, `in_start`, `out_list` and `in_list`, or calls [[ns_finish_lists]] for them.
subroutine ns_allocate(s, n, m)
    type(ns_t), intent(out) :: s    !! the graph
    integer, intent(in) :: n        !! number of nodes
    integer, intent(in) :: m        !! number of edges

    s%n = n
    s%m = m
    allocate (s%tail(m), s%head(m), s%minlen(m), s%weight(m), source=0)
    allocate (s%rank(n), source=0)
    allocate (s%normal(n), source=.true.)
    allocate (s%out_start(n + 1), s%in_start(n + 1), source=1)
    allocate (s%out_list(m), s%in_list(m), source=0)
end subroutine ns_allocate

!> Build the lists of edges of each node from `tail` and `head`, in the order of the edges.
subroutine ns_finish_lists(s)
    type(ns_t), intent(inout) :: s  !! the graph
    integer, allocatable :: out_next(:), in_next(:)
    integer :: e, v

    s%out_start = 0
    s%in_start = 0
    do e = 1, s%m
        s%out_start(s%tail(e)) = s%out_start(s%tail(e)) + 1
        s%in_start(s%head(e)) = s%in_start(s%head(e)) + 1
    end do
    allocate (out_next(s%n + 1), in_next(s%n + 1))
    out_next(1) = 1
    in_next(1) = 1
    do v = 1, s%n
        out_next(v + 1) = out_next(v) + s%out_start(v)
        in_next(v + 1) = in_next(v) + s%in_start(v)
    end do
    s%out_start = out_next
    s%in_start = in_next
    do e = 1, s%m
        s%out_list(out_next(s%tail(e))) = e
        out_next(s%tail(e)) = out_next(s%tail(e)) + 1
        s%in_list(in_next(s%head(e))) = e
        in_next(s%head(e)) = in_next(s%head(e)) + 1
    end do
end subroutine ns_finish_lists

!> Solve. `stat` is 0, or 1 if the graph is not connected, or 2 if something is seriously wrong.
subroutine ns_solve(s, balance, maxiter, stat)
    type(ns_t), intent(inout) :: s      !! the graph, with its result in `rank`
    integer, intent(in) :: balance      !! 1 to move nodes to less populated ranks, 2 to centre, or 0
    integer, intent(in) :: maxiter      !! the limit of iterations
    integer, intent(out) :: stat        !! the result
    logical :: feasible
    integer :: e, f, iter

    stat = 0
    if (s%n == 0) return
    call init_graph(s, feasible)
    if (.not. feasible) call init_rank(s, stat)
    if (stat /= 0) return
    call feasible_tree(s, stat)
    if (stat /= 0) return
    if (maxiter <= 0) return

    iter = 0
    do
        e = leave_edge(s)
        if (e == 0) exit
        f = enter_edge(s, e)
        if (f == 0) then
            stat = 2
            return
        end if
        call update(s, e, f, stat)
        if (stat /= 0) return
        iter = iter + 1
        if (iter >= maxiter) exit
    end do
    select case (balance)
    case (1)
        call tb_balance(s)
    case (2)
        call lr_balance(s)
    case default
        block
            integer :: maxrank
            call scan_and_normalize(s, maxrank)
        end block
    end select
end subroutine ns_solve

!> The length of edge `e` beyond its smallest.
integer function slack(s, e)
    type(ns_t), intent(in) :: s     !! the graph
    integer, intent(in) :: e        !! the edge

    slack = s%rank(s%head(e)) - s%rank(s%tail(e)) - s%minlen(e)
end function slack

!> Whether `a <= b <= c`.
logical function seq(a, b, c)
    integer, intent(in) :: a    !! the lower limit
    integer, intent(in) :: b    !! the value
    integer, intent(in) :: c    !! the upper limit

    seq = a <= b .and. b <= c
end function seq

!> Whether `e` is in the tree.
logical function in_tree(s, e)
    type(ns_t), intent(in) :: s     !! the graph
    integer, intent(in) :: e        !! the edge

    in_tree = s%tree_index(e) > 0
end function in_tree

!> Set up what the solver keeps, and report whether the ranks given are feasible.
subroutine init_graph(s, feasible)
    type(ns_t), intent(inout) :: s  !! the graph
    logical, intent(out) :: feasible    !! the result
    integer :: e, v

    if (allocated(s%cutvalue)) deallocate (s%cutvalue)
    allocate (s%cutvalue(s%m), source=0_int64)
    call reset(s%tree_index, s%m)
    call reset(s%tree_edge, s%n)
    call reset(s%tree_out, s%m)
    call reset(s%tree_in, s%m)
    call reset(s%tree_out_size, s%n)
    call reset(s%tree_in_size, s%n)
    call reset(s%par, s%n)
    call reset(s%low, s%n)
    call reset(s%lim, s%n)
    call reset(s%priority, s%n)
    call reset(s%subtree, s%n)
    call reset(s%st_rep, s%n)
    call reset(s%st_size, s%n)
    call reset(s%st_heap, s%n)
    call reset(s%st_par, s%n)
    if (allocated(s%heap)) deallocate (s%heap)
    allocate (s%heap(0:s%n), source=0)
    if (allocated(s%stack)) deallocate (s%stack)
    allocate (s%stack(5, s%n + 1), source=0)
    s%tree_count = 0
    s%search_index = 0
    s%search_size = SEARCHSIZE

    feasible = .true.
    do v = 1, s%n
        s%priority(v) = s%in_start(v + 1) - s%in_start(v)
    end do
    do e = 1, s%m
        if (s%rank(s%head(e)) - s%rank(s%tail(e)) < s%minlen(e)) feasible = .false.
    end do
end subroutine init_graph

!> Make `array` a zeroed array of `n` values.
subroutine reset(array, n)
    integer, allocatable, intent(inout) :: array(:) !! the array
    integer, intent(in) :: n                        !! its size

    if (allocated(array)) deallocate (array)
    allocate (array(n), source=0)
end subroutine reset

!> Add `e` to the tree.
subroutine add_tree_edge(s, e)
    type(ns_t), intent(inout) :: s  !! the graph
    integer, intent(in) :: e        !! the edge
    integer :: v

    s%tree_count = s%tree_count + 1
    s%tree_index(e) = s%tree_count
    s%tree_edge(s%tree_count) = e
    v = s%tail(e)
    s%tree_out(s%out_start(v) + s%tree_out_size(v)) = e
    s%tree_out_size(v) = s%tree_out_size(v) + 1
    v = s%head(e)
    s%tree_in(s%in_start(v) + s%tree_in_size(v)) = e
    s%tree_in_size(v) = s%tree_in_size(v) + 1
end subroutine add_tree_edge

!> Replace tree edge `e` by `f`.
subroutine exchange_tree_edges(s, e, f)
    type(ns_t), intent(inout) :: s  !! the graph
    integer, intent(in) :: e        !! the edge that leaves
    integer, intent(in) :: f        !! the edge that enters
    integer :: v, j, last

    s%tree_index(f) = s%tree_index(e)
    s%tree_edge(s%tree_index(e)) = f
    s%tree_index(e) = 0

    v = s%tail(e)
    s%tree_out_size(v) = s%tree_out_size(v) - 1
    last = s%out_start(v) + s%tree_out_size(v)
    do j = s%out_start(v), last
        if (s%tree_out(j) == e) exit
    end do
    s%tree_out(j) = s%tree_out(last)
    s%tree_out(last) = 0
    v = s%head(e)
    s%tree_in_size(v) = s%tree_in_size(v) - 1
    last = s%in_start(v) + s%tree_in_size(v)
    do j = s%in_start(v), last
        if (s%tree_in(j) == e) exit
    end do
    s%tree_in(j) = s%tree_in(last)
    s%tree_in(last) = 0

    v = s%tail(f)
    s%tree_out(s%out_start(v) + s%tree_out_size(v)) = f
    s%tree_out_size(v) = s%tree_out_size(v) + 1
    v = s%head(f)
    s%tree_in(s%in_start(v) + s%tree_in_size(v)) = f
    s%tree_in_size(v) = s%tree_in_size(v) + 1
end subroutine exchange_tree_edges

!> Invalidate what the depth-first search keeps, walking up the tree from `to_node` to `lca`. This
!> improves the pruning in [[dfs_range]].
subroutine invalidate_path(s, lca, to_node)
    type(ns_t), intent(inout) :: s      !! the graph
    integer, intent(in) :: lca          !! the lowest common ancestor
    integer, intent(in) :: to_node      !! where to start
    integer :: v, e

    v = to_node
    do
        if (s%low(v) == -1) exit
        s%low(v) = -1
        e = s%par(v)
        if (e == 0) exit
        if (s%lim(v) >= s%lim(lca)) exit
        if (s%lim(s%tail(e)) > s%lim(s%head(e))) then
            v = s%tail(e)
        else
            v = s%head(e)
        end if
    end do
end subroutine invalidate_path

!> A first feasible ranking, from a topological sort.
subroutine init_rank(s, stat)
    type(ns_t), intent(inout) :: s      !! the graph
    integer, intent(out) :: stat        !! 2 if the graph has a cycle
    integer, allocatable :: q(:)
    integer :: v, i, e, first, last, ctr

    allocate (q(s%n))
    first = 1
    last = 0
    ctr = 0
    do v = 1, s%n
        if (s%priority(v) == 0) then
            last = last + 1
            q(last) = v
        end if
    end do
    do while (first <= last)
        v = q(first)
        first = first + 1
        s%rank(v) = 0
        ctr = ctr + 1
        do i = s%in_start(v), s%in_start(v + 1) - 1
            e = s%in_list(i)
            s%rank(v) = max(s%rank(v), s%rank(s%tail(e)) + s%minlen(e))
        end do
        do i = s%out_start(v), s%out_start(v + 1) - 1
            e = s%out_list(i)
            s%priority(s%head(e)) = s%priority(s%head(e)) - 1
            if (s%priority(s%head(e)) == 0) then
                last = last + 1
                q(last) = s%head(e)
            end if
        end do
    end do
    stat = 0
    if (ctr /= s%n) stat = 2
end subroutine init_rank

!> The tree edge to remove, with a negative cut value, or 0 if there is none.
integer function leave_edge(s) result(rv)
    type(ns_t), intent(inout) :: s  !! the graph
    integer :: f, j, cnt

    rv = 0
    cnt = 0
    j = s%search_index
    do while (s%search_index < s%tree_count)
        f = s%tree_edge(s%search_index + 1)
        if (s%cutvalue(f) < 0) then
            if (rv /= 0) then
                if (s%cutvalue(rv) > s%cutvalue(f)) rv = f
            else
                rv = f
            end if
            cnt = cnt + 1
            if (cnt >= s%search_size) return
        end if
        s%search_index = s%search_index + 1
    end do
    if (j > 0) then
        s%search_index = 0
        do while (s%search_index < j)
            f = s%tree_edge(s%search_index + 1)
            if (s%cutvalue(f) < 0) then
                if (rv /= 0) then
                    if (s%cutvalue(rv) > s%cutvalue(f)) rv = f
                else
                    rv = f
                end if
                cnt = cnt + 1
                if (cnt >= s%search_size) return
            end if
            s%search_index = s%search_index + 1
        end do
    end if
end function leave_edge

!> The edge to enter the tree in place of `e`.
integer function enter_edge(s, e)
    type(ns_t), intent(inout) :: s  !! the graph
    integer, intent(in) :: e        !! the edge that leaves

    if (s%lim(s%tail(e)) < s%lim(s%head(e))) then
        enter_edge = dfs_enter_inedge(s, s%tail(e), s%low(s%tail(e)), s%lim(s%tail(e)))
    else
        enter_edge = dfs_enter_outedge(s, s%head(e), s%low(s%head(e)), s%lim(s%head(e)))
    end if
end function enter_edge

!> The tightest edge leaving the subtree of `start`.
integer function dfs_enter_outedge(s, start, low, lim) result(enter)
    type(ns_t), intent(inout) :: s  !! the graph
    integer, intent(in) :: start    !! the root of the subtree
    integer, intent(in) :: low      !! the lowest search index of the subtree
    integer, intent(in) :: lim      !! the highest search index of the subtree
    integer :: v, e, i, sp, best, this

    enter = 0
    best = huge(best)
    sp = 1
    s%stack(1, 1) = start
    do while (sp > 0)
        v = s%stack(1, sp)
        sp = sp - 1
        do i = s%out_start(v), s%out_start(v + 1) - 1
            e = s%out_list(i)
            if (.not. in_tree(s, e)) then
                if (.not. seq(low, s%lim(s%head(e)), lim)) then
                    this = slack(s, e)
                    if (this < best .or. enter == 0) then
                        enter = e
                        best = this
                    end if
                end if
            else if (s%lim(s%head(e)) < s%lim(v)) then
                sp = sp + 1
                s%stack(1, sp) = s%head(e)
            end if
        end do
        do i = s%in_start(v), s%in_start(v) + s%tree_in_size(v) - 1
            if (best <= 0) exit
            e = s%tree_in(i)
            if (s%lim(s%tail(e)) < s%lim(v)) then
                sp = sp + 1
                s%stack(1, sp) = s%tail(e)
            end if
        end do
    end do
end function dfs_enter_outedge

!> The tightest edge entering the subtree of `start`.
integer function dfs_enter_inedge(s, start, low, lim) result(enter)
    type(ns_t), intent(inout) :: s  !! the graph
    integer, intent(in) :: start    !! the root of the subtree
    integer, intent(in) :: low      !! the lowest search index of the subtree
    integer, intent(in) :: lim      !! the highest search index of the subtree
    integer :: v, e, i, sp, best, this

    enter = 0
    best = huge(best)
    sp = 1
    s%stack(1, 1) = start
    do while (sp > 0)
        v = s%stack(1, sp)
        sp = sp - 1
        do i = s%in_start(v), s%in_start(v + 1) - 1
            e = s%in_list(i)
            if (.not. in_tree(s, e)) then
                if (.not. seq(low, s%lim(s%tail(e)), lim)) then
                    this = slack(s, e)
                    if (this < best .or. enter == 0) then
                        enter = e
                        best = this
                    end if
                end if
            else if (s%lim(s%tail(e)) < s%lim(v)) then
                sp = sp + 1
                s%stack(1, sp) = s%tail(e)
            end if
        end do
        do i = s%out_start(v), s%out_start(v) + s%tree_out_size(v) - 1
            if (best <= 0) exit
            e = s%tree_out(i)
            if (s%lim(s%head(e)) < s%lim(v)) then
                sp = sp + 1
                s%stack(1, sp) = s%head(e)
            end if
        end do
    end do
end function dfs_enter_inedge

!> Find the initial tight subtree holding `start`, and return how many nodes it has.
subroutine tight_subtree_search(s, start, st, rv)
    type(ns_t), intent(inout) :: s      !! the graph
    integer, intent(in) :: start        !! where to start
    integer, intent(in) :: st           !! the subtree
    integer, intent(out) :: rv          !! the size of the subtree
    integer :: sp, v, e, last

    ! stack: 1 is the node, 2 the position among in-edges, 3 among out-edges, 4 the result
    s%subtree(start) = st
    sp = 1
    s%stack(:, 1) = [start, 0, 0, 1, 0]
    outer: do while (sp > 0)
        v = s%stack(1, sp)
        do while (s%in_start(v) + s%stack(2, sp) < s%in_start(v + 1))
            e = s%in_list(s%in_start(v) + s%stack(2, sp))
            s%stack(2, sp) = s%stack(2, sp) + 1
            if (in_tree(s, e)) cycle
            if (s%subtree(s%tail(e)) == 0 .and. slack(s, e) == 0) then
                call add_tree_edge(s, e)
                s%subtree(s%tail(e)) = st
                sp = sp + 1
                s%stack(:, sp) = [s%tail(e), 0, 0, 1, 0]
                cycle outer
            end if
        end do
        do while (s%out_start(v) + s%stack(3, sp) < s%out_start(v + 1))
            e = s%out_list(s%out_start(v) + s%stack(3, sp))
            s%stack(3, sp) = s%stack(3, sp) + 1
            if (in_tree(s, e)) cycle
            if (s%subtree(s%head(e)) == 0 .and. slack(s, e) == 0) then
                call add_tree_edge(s, e)
                s%subtree(s%head(e)) = st
                sp = sp + 1
                s%stack(:, sp) = [s%head(e), 0, 0, 1, 0]
                cycle outer
            end if
        end do
        last = s%stack(4, sp)
        sp = sp - 1
        if (sp == 0) then
            rv = last
        else
            s%stack(4, sp) = s%stack(4, sp) + last
        end if
    end do outer
end subroutine tight_subtree_search

!> The root of the set of subtrees holding node `n0`.
integer function st_find(s, n0) result(s0)
    type(ns_t), intent(inout) :: s  !! the graph
    integer, intent(in) :: n0       !! the node

    s0 = s%subtree(n0)
    do while (s%st_par(s0) /= 0 .and. s%st_par(s0) /= s0)
        if (s%st_par(s%st_par(s0)) /= 0) s%st_par(s0) = s%st_par(s%st_par(s0))
        s0 = s%st_par(s0)
    end do
end function st_find

!> Join two sets of subtrees, and return the root.
integer function st_union(s, s0, s1) result(r)
    type(ns_t), intent(inout) :: s  !! the graph
    integer, intent(in) :: s0       !! one set
    integer, intent(in) :: s1       !! the other set
    integer :: r0, r1

    r0 = s0
    do while (s%st_par(r0) /= 0 .and. s%st_par(r0) /= r0)
        r0 = s%st_par(r0)
    end do
    r1 = s1
    do while (s%st_par(r1) /= 0 .and. s%st_par(r1) /= r1)
        r1 = s%st_par(r1)
    end do
    if (r0 == r1) then
        r = r0
        return
    end if
    if (s%st_heap(r1) < 0) then
        r = r0
    else if (s%st_heap(r0) < 0) then
        r = r1
    else if (s%st_size(r1) < s%st_size(r0)) then
        r = r0
    else
        r = r1
    end if
    s%st_par(r0) = r
    s%st_par(r1) = r
    s%st_size(r) = s%st_size(r0) + s%st_size(r1)
end function st_union

!> The tightest edge from the tree of `start` to another tree.
integer function inter_tree_edge(s, start) result(best)
    type(ns_t), intent(inout) :: s  !! the graph
    integer, intent(in) :: start    !! a node of the tree
    integer :: sp, v, from, e, ts

    ! stack: 1 is the node, 2 its tree, 3 where the search came from, 4 and 5 the positions
    best = 0
    sp = 1
    s%stack(:, 1) = [start, st_find(s, start), 0, 0, 0]
    outer: do while (sp > 0)
        if (s%stack(4, sp) == 0 .and. s%stack(5, sp) == 0 .and. best /= 0) then
            if (slack(s, best) == 0) then
                sp = sp - 1
                cycle
            end if
        end if
        v = s%stack(1, sp)
        ts = s%stack(2, sp)
        from = s%stack(3, sp)
        do while (s%out_start(v) + s%stack(4, sp) < s%out_start(v + 1))
            e = s%out_list(s%out_start(v) + s%stack(4, sp))
            s%stack(4, sp) = s%stack(4, sp) + 1
            if (in_tree(s, e)) then
                if (s%head(e) == from) cycle
                sp = sp + 1
                s%stack(:, sp) = [s%head(e), st_find(s, s%head(e)), v, 0, 0]
                cycle outer
            else if (st_find(s, s%head(e)) /= ts) then
                if (best == 0) then
                    best = e
                else if (slack(s, e) < slack(s, best)) then
                    best = e
                end if
            end if
        end do
        do while (s%in_start(v) + s%stack(5, sp) < s%in_start(v + 1))
            e = s%in_list(s%in_start(v) + s%stack(5, sp))
            s%stack(5, sp) = s%stack(5, sp) + 1
            if (in_tree(s, e)) then
                if (s%tail(e) == from) cycle
                sp = sp + 1
                s%stack(:, sp) = [s%tail(e), st_find(s, s%tail(e)), v, 0, 0]
                cycle outer
            else if (st_find(s, s%tail(e)) /= ts) then
                if (best == 0) then
                    best = e
                else if (slack(s, e) < slack(s, best)) then
                    best = e
                end if
            end if
        end do
        sp = sp - 1
    end do outer
end function inter_tree_edge

!> Restore the heap order below position `start`.
subroutine st_heapify(s, start)
    type(ns_t), intent(inout) :: s  !! the graph
    integer, intent(in) :: start    !! the position
    integer :: i, left, right, smallest, t

    i = start
    do
        left = 2*(i + 1) - 1
        right = 2*(i + 1)
        smallest = i
        if (left < s%heap_size) then
            if (s%st_size(s%heap(left)) < s%st_size(s%heap(smallest))) smallest = left
        end if
        if (right < s%heap_size) then
            if (s%st_size(s%heap(right)) < s%st_size(s%heap(smallest))) smallest = right
        end if
        if (smallest == i) exit
        t = s%heap(i)
        s%heap(i) = s%heap(smallest)
        s%heap(smallest) = t
        s%st_heap(s%heap(i)) = i
        s%st_heap(s%heap(smallest)) = smallest
        i = smallest
        if (i >= s%heap_size) exit
    end do
end subroutine st_heapify

!> Move the tree around `start` by `delta` ranks.
subroutine tree_adjust(s, start, delta)
    type(ns_t), intent(inout) :: s  !! the graph
    integer, intent(in) :: start    !! a node of the tree
    integer, intent(in) :: delta    !! how far
    integer :: sp, v, from, e, w

    ! stack: 1 is the node, 2 where the search came from, 3 and 4 the positions
    sp = 1
    s%stack(:, 1) = [start, 0, 0, 0, 0]
    outer: do while (sp > 0)
        v = s%stack(1, sp)
        from = s%stack(2, sp)
        if (s%stack(3, sp) == 0 .and. s%stack(4, sp) == 0) s%rank(v) = s%rank(v) + delta
        do while (s%stack(3, sp) < s%tree_in_size(v))
            e = s%tree_in(s%in_start(v) + s%stack(3, sp))
            s%stack(3, sp) = s%stack(3, sp) + 1
            w = s%tail(e)
            if (w /= from) then
                sp = sp + 1
                s%stack(:, sp) = [w, v, 0, 0, 0]
                cycle outer
            end if
        end do
        do while (s%stack(4, sp) < s%tree_out_size(v))
            e = s%tree_out(s%out_start(v) + s%stack(4, sp))
            s%stack(4, sp) = s%stack(4, sp) + 1
            w = s%head(e)
            if (w /= from) then
                sp = sp + 1
                s%stack(:, sp) = [w, v, 0, 0, 0]
                cycle outer
            end if
        end do
        sp = sp - 1
    end do outer
end subroutine tree_adjust

!> Join the two trees at the ends of `e`, and return the result.
integer function merge_trees(s, e)
    type(ns_t), intent(inout) :: s  !! the graph
    integer, intent(in) :: e        !! the edge that enters the tree
    integer :: t0, t1, delta

    t0 = st_find(s, s%tail(e))
    t1 = st_find(s, s%head(e))
    if (s%st_heap(t0) < 0) then
        delta = slack(s, e)
        if (delta /= 0) call tree_adjust(s, s%st_rep(t0), delta)
    else
        delta = -slack(s, e)
        if (delta /= 0) call tree_adjust(s, s%st_rep(t1), delta)
    end if
    call add_tree_edge(s, e)
    merge_trees = st_union(s, t0, t1)
end function merge_trees

!> Build the initial tight tree, moving ranks as needed. `stat` is 1 if the graph is not connected.
subroutine feasible_tree(s, stat)
    type(ns_t), intent(inout) :: s      !! the graph
    integer, intent(out) :: stat        !! the result
    integer :: v, count_trees, i, tree0, tree1, ee, size

    stat = 0
    s%subtree = 0
    count_trees = 0
    do v = 1, s%n
        if (s%subtree(v) == 0) then
            count_trees = count_trees + 1
            s%st_rep(count_trees) = v
            call tight_subtree_search(s, v, count_trees, size)
            s%st_size(count_trees) = size
            s%st_par(count_trees) = count_trees
        end if
    end do

    s%heap_size = count_trees
    do i = 0, count_trees - 1
        s%heap(i) = i + 1
        s%st_heap(i + 1) = i
    end do
    do i = s%heap_size/2, 0, -1
        call st_heapify(s, i)
    end do
    do while (s%heap_size > 1)
        tree0 = s%heap(0)
        s%st_heap(tree0) = -1
        s%heap(0) = s%heap(s%heap_size - 1)
        s%st_heap(s%heap(0)) = 0
        s%heap(s%heap_size - 1) = tree0
        s%heap_size = s%heap_size - 1
        call st_heapify(s, 0)
        ee = inter_tree_edge(s, s%st_rep(tree0))
        if (ee == 0) then
            stat = 1
            return
        end if
        tree1 = merge_trees(s, ee)
        call st_heapify(s, s%st_heap(tree1))
    end do
    if (s%tree_count /= s%n - 1) then
        stat = 2
        return
    end if
    call dfs_range_init(s, 1)
    call dfs_cutval(s, 1)
end subroutine feasible_tree

!> Walk up from `start` to the common ancestor with `w`, setting new cut values.
integer function treeupdate(s, start, w, cutvalue, dir) result(v)
    type(ns_t), intent(inout) :: s          !! the graph
    integer, intent(in) :: start            !! where to start
    integer, intent(in) :: w                !! the other end
    integer(int64), intent(in) :: cutvalue  !! the change
    logical, intent(in) :: dir              !! its direction
    integer :: e
    logical :: d

    v = start
    do while (.not. seq(s%low(v), s%lim(w), s%lim(v)))
        e = s%par(v)
        if (v == s%tail(e)) then
            d = dir
        else
            d = .not. dir
        end if
        if (d) then
            s%cutvalue(e) = s%cutvalue(e) + cutvalue
        else
            s%cutvalue(e) = s%cutvalue(e) - cutvalue
        end if
        if (s%lim(s%tail(e)) > s%lim(s%head(e))) then
            v = s%tail(e)
        else
            v = s%head(e)
        end if
    end do
end function treeupdate

!> Move the part of the tree below `start` by `delta` ranks.
subroutine rerank(s, start, delta)
    type(ns_t), intent(inout) :: s  !! the graph
    integer, intent(in) :: start    !! where to start
    integer, intent(in) :: delta    !! how far
    integer :: sp, v, e

    ! stack: 1 is the node, 2 and 3 the positions
    sp = 1
    s%stack(:, 1) = [start, 0, 0, 0, 0]
    outer: do while (sp > 0)
        v = s%stack(1, sp)
        if (s%stack(2, sp) == 0 .and. s%stack(3, sp) == 0) s%rank(v) = s%rank(v) - delta
        do while (s%stack(2, sp) < s%tree_out_size(v))
            e = s%tree_out(s%out_start(v) + s%stack(2, sp))
            s%stack(2, sp) = s%stack(2, sp) + 1
            if (e /= s%par(v)) then
                sp = sp + 1
                s%stack(:, sp) = [s%head(e), 0, 0, 0, 0]
                cycle outer
            end if
        end do
        do while (s%stack(3, sp) < s%tree_in_size(v))
            e = s%tree_in(s%in_start(v) + s%stack(3, sp))
            s%stack(3, sp) = s%stack(3, sp) + 1
            if (e /= s%par(v)) then
                sp = sp + 1
                s%stack(:, sp) = [s%tail(e), 0, 0, 0, 0]
                cycle outer
            end if
        end do
        sp = sp - 1
    end do outer
end subroutine rerank

!> Exchange tree edge `e` with non-tree edge `f`, and update ranks and cut values.
subroutine update(s, e, f, stat)
    type(ns_t), intent(inout) :: s      !! the graph
    integer, intent(in) :: e            !! the edge that leaves
    integer, intent(in) :: f            !! the edge that enters
    integer, intent(out) :: stat        !! 2 if the tree is inconsistent
    integer(int64) :: cutvalue
    integer :: delta, count, lca, lca_low, other

    stat = 0
    delta = slack(s, f)
    if (delta > 0) then
        count = s%tree_in_size(s%tail(e)) + s%tree_out_size(s%tail(e))
        if (count == 1) then
            call rerank(s, s%tail(e), delta)
        else
            count = s%tree_in_size(s%head(e)) + s%tree_out_size(s%head(e))
            if (count == 1) then
                call rerank(s, s%head(e), -delta)
            else if (s%lim(s%tail(e)) < s%lim(s%head(e))) then
                call rerank(s, s%tail(e), delta)
            else
                call rerank(s, s%head(e), -delta)
            end if
        end if
    end if

    cutvalue = s%cutvalue(e)
    lca = treeupdate(s, s%tail(f), s%head(f), cutvalue, .true.)
    other = treeupdate(s, s%head(f), s%tail(f), cutvalue, .false.)
    if (other /= lca) then
        stat = 2
        return
    end if

    lca_low = s%low(lca)
    call invalidate_path(s, lca, s%head(f))
    call invalidate_path(s, lca, s%tail(f))

    s%cutvalue(f) = -cutvalue
    s%cutvalue(e) = 0
    call exchange_tree_edges(s, e, f)
    call dfs_range(s, lca, s%par(lca), lca_low)
end subroutine update

!> Make the smallest rank of the nodes of the drawing 0, and return the largest.
subroutine scan_and_normalize(s, maxrank)
    type(ns_t), intent(inout) :: s      !! the graph
    integer, intent(out) :: maxrank     !! the largest rank
    integer :: minrank, v

    minrank = huge(minrank)
    maxrank = -huge(maxrank)
    do v = 1, s%n
        if (s%normal(v)) then
            minrank = min(minrank, s%rank(v))
            maxrank = max(maxrank, s%rank(v))
        end if
    end do
    s%rank = s%rank - minrank
    maxrank = maxrank - minrank
end subroutine scan_and_normalize

!> Centre what can move, for x coordinates.
subroutine lr_balance(s)
    type(ns_t), intent(inout) :: s  !! the graph
    integer :: i, e, f, delta

    do i = 1, s%tree_count
        e = s%tree_edge(i)
        if (s%cutvalue(e) == 0) then
            f = enter_edge(s, e)
            if (f == 0) cycle
            delta = slack(s, f)
            if (delta <= 1) cycle
            if (s%lim(s%tail(e)) < s%lim(s%head(e))) then
                call rerank(s, s%tail(e), delta/2)
            else
                call rerank(s, s%head(e), -delta/2)
            end if
        end if
    end do
end subroutine lr_balance

!> Move nodes that are not tight to less populated ranks.
subroutine tb_balance(s)
    type(ns_t), intent(inout) :: s  !! the graph
    integer, allocatable :: nrank(:), order(:), bucket(:)
    integer :: maxrank, v, i, ii, e, inweight, outweight, low, high, choice, lowest, r

    call scan_and_normalize(s, maxrank)
    allocate (nrank(0:maxrank), source=0)

    ! the nodes by increasing rank, in a stable order, as the same input must give the same output
    ! on every platform, which the `qsort` of the original does not promise
    allocate (order(s%n))
    lowest = minval(s%rank)
    allocate (bucket(lowest:maxval(s%rank) + 1), source=0)
    do v = 1, s%n
        bucket(s%rank(v)) = bucket(s%rank(v)) + 1
    end do
    i = 1
    do r = lowest, ubound(bucket, 1)
        ii = bucket(r)
        bucket(r) = i
        i = i + ii
    end do
    do v = 1, s%n
        order(bucket(s%rank(v))) = v
        bucket(s%rank(v)) = bucket(s%rank(v)) + 1
    end do

    do ii = 1, s%n
        v = order(ii)
        if (s%normal(v)) nrank(s%rank(v)) = nrank(s%rank(v)) + 1
    end do
    do ii = 1, s%n
        v = order(ii)
        if (.not. s%normal(v)) cycle
        inweight = 0
        outweight = 0
        low = 0
        high = maxrank
        do i = s%in_start(v), s%in_start(v + 1) - 1
            e = s%in_list(i)
            inweight = inweight + s%weight(e)
            low = max(low, s%rank(s%tail(e)) + s%minlen(e))
        end do
        do i = s%out_start(v), s%out_start(v + 1) - 1
            e = s%out_list(i)
            outweight = outweight + s%weight(e)
            high = min(high, s%rank(s%head(e)) - s%minlen(e))
        end do
        if (low < 0) low = 0
        if (inweight == outweight) then
            choice = low
            do i = low + 1, high
                if (nrank(i) < nrank(choice)) choice = i
            end do
            nrank(s%rank(v)) = nrank(s%rank(v)) - 1
            nrank(choice) = nrank(choice) + 1
            s%rank(v) = choice
        end if
    end do
end subroutine tb_balance

!> Set the cut value of `f`, assuming the values of the edges on one side are set.
subroutine x_cutval(s, f)
    type(ns_t), intent(inout) :: s  !! the graph
    integer, intent(in) :: f        !! the edge
    integer(int64) :: sum
    integer :: v, dir, i

    if (s%par(s%tail(f)) == f) then
        v = s%tail(f)
        dir = 1
    else
        v = s%head(f)
        dir = -1
    end if
    sum = 0
    do i = s%out_start(v), s%out_start(v + 1) - 1
        sum = sum + x_val(s, s%out_list(i), v, dir)
    end do
    do i = s%in_start(v), s%in_start(v + 1) - 1
        sum = sum + x_val(s, s%in_list(i), v, dir)
    end do
    s%cutvalue(f) = sum
end subroutine x_cutval

!> What edge `e` of node `v` adds to a cut value.
integer(int64) function x_val(s, e, v, dir) result(rv)
    type(ns_t), intent(in) :: s     !! the graph
    integer, intent(in) :: e        !! the edge
    integer, intent(in) :: v        !! the node
    integer, intent(in) :: dir      !! which side has been searched
    integer :: other, d
    logical :: f

    if (s%tail(e) == v) then
        other = s%head(e)
    else
        other = s%tail(e)
    end if
    if (.not. seq(s%low(v), s%lim(other), s%lim(v))) then
        f = .true.
        rv = s%weight(e)
    else
        f = .false.
        if (in_tree(s, e)) then
            rv = s%cutvalue(e)
        else
            rv = 0
        end if
        rv = rv - s%weight(e)
    end if
    if (dir > 0) then
        d = merge(1, -1, s%head(e) == v)
    else
        d = merge(1, -1, s%tail(e) == v)
    end if
    if (f) d = -d
    if (d < 0) rv = -rv
end function x_val

!> Set all cut values, from the leaves of the tree up.
subroutine dfs_cutval(s, start)
    type(ns_t), intent(inout) :: s  !! the graph
    integer, intent(in) :: start    !! the root
    integer :: sp, v, par, e

    ! stack: 1 is the node, 2 its parent edge, 3 and 4 the positions
    sp = 1
    s%stack(:, 1) = [start, 0, 0, 0, 0]
    outer: do while (sp > 0)
        v = s%stack(1, sp)
        par = s%stack(2, sp)
        do while (s%stack(3, sp) < s%tree_out_size(v))
            e = s%tree_out(s%out_start(v) + s%stack(3, sp))
            s%stack(3, sp) = s%stack(3, sp) + 1
            if (e /= par) then
                sp = sp + 1
                s%stack(:, sp) = [s%head(e), e, 0, 0, 0]
                cycle outer
            end if
        end do
        do while (s%stack(4, sp) < s%tree_in_size(v))
            e = s%tree_in(s%in_start(v) + s%stack(4, sp))
            s%stack(4, sp) = s%stack(4, sp) + 1
            if (e /= par) then
                sp = sp + 1
                s%stack(:, sp) = [s%tail(e), e, 0, 0, 0]
                cycle outer
            end if
        end do
        if (par /= 0) call x_cutval(s, par)
        sp = sp - 1
    end do outer
end subroutine dfs_cutval

!> Set the parent edge and search range of every node.
subroutine dfs_range_init(s, start)
    type(ns_t), intent(inout) :: s  !! the graph
    integer, intent(in) :: start    !! the root
    integer :: sp, v, par, e, n, lim

    ! stack: 1 is the node, 2 its parent edge, 3 the range, 4 and 5 the positions
    s%par(start) = 0
    s%low(start) = 1
    sp = 1
    s%stack(:, 1) = [start, 0, 1, 0, 0]
    outer: do while (sp > 0)
        v = s%stack(1, sp)
        par = s%stack(2, sp)
        do while (s%stack(4, sp) < s%tree_out_size(v))
            e = s%tree_out(s%out_start(v) + s%stack(4, sp))
            s%stack(4, sp) = s%stack(4, sp) + 1
            if (e /= par) then
                n = s%head(e)
                s%par(n) = e
                s%low(n) = s%stack(3, sp)
                s%stack(:, sp + 1) = [n, e, s%stack(3, sp), 0, 0]
                sp = sp + 1
                cycle outer
            end if
        end do
        do while (s%stack(5, sp) < s%tree_in_size(v))
            e = s%tree_in(s%in_start(v) + s%stack(5, sp))
            s%stack(5, sp) = s%stack(5, sp) + 1
            if (e /= par) then
                n = s%tail(e)
                s%par(n) = e
                s%low(n) = s%stack(3, sp)
                s%stack(:, sp + 1) = [n, e, s%stack(3, sp), 0, 0]
                sp = sp + 1
                cycle outer
            end if
        end do
        lim = s%stack(3, sp)
        s%lim(v) = lim
        sp = sp - 1
        if (sp > 0) s%stack(3, sp) = lim + 1
    end do outer
end subroutine dfs_range_init

!> Update the parent edges and search ranges below `start`, stopping where nothing has changed.
subroutine dfs_range(s, start, start_par, start_low)
    type(ns_t), intent(inout) :: s      !! the graph
    integer, intent(in) :: start        !! where to start
    integer, intent(in) :: start_par    !! its parent edge
    integer, intent(in) :: start_low    !! its lowest search index
    integer :: sp, v, par, e, n, lim

    if (s%par(start) == start_par .and. s%low(start) == start_low) return
    s%par(start) = start_par
    s%low(start) = start_low
    sp = 1
    s%stack(:, 1) = [start, start_par, start_low, 0, 0]
    outer: do while (sp > 0)
        v = s%stack(1, sp)
        par = s%stack(2, sp)
        do while (s%stack(4, sp) < s%tree_out_size(v))
            e = s%tree_out(s%out_start(v) + s%stack(4, sp))
            s%stack(4, sp) = s%stack(4, sp) + 1
            if (e /= par) then
                n = s%head(e)
                call visit(n, e)
                cycle outer
            end if
        end do
        do while (s%stack(5, sp) < s%tree_in_size(v))
            e = s%tree_in(s%in_start(v) + s%stack(5, sp))
            s%stack(5, sp) = s%stack(5, sp) + 1
            if (e /= par) then
                n = s%tail(e)
                call visit(n, e)
                cycle outer
            end if
        end do
        lim = s%stack(3, sp)
        s%lim(v) = lim
        sp = sp - 1
        if (sp > 0) s%stack(3, sp) = lim + 1
    end do outer

contains

    !> Go down to child `child`, or step over it if it has not changed.
    subroutine visit(child, edge)
        integer, intent(in) :: child    !! the child
        integer, intent(in) :: edge     !! the edge to it

        if (s%par(child) == edge .and. s%low(child) == s%stack(3, sp)) then
            s%stack(3, sp) = s%lim(child) + 1
        else
            s%par(child) = edge
            s%low(child) = s%stack(3, sp)
            s%stack(:, sp + 1) = [child, edge, s%stack(3, sp), 0, 0]
            sp = sp + 1
        end if
    end subroutine visit
end subroutine dfs_range

end module fdot_simplex
