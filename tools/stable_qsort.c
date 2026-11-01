/* A stable merge sort in place of `qsort`, for comparing with `dot`.
 *
 * `dot` leaves some ties to `qsort`, in `TB_balance` of ns.c and `edgecmp` of dotsplines.c, so its result
 * depends on the C library. fdot sorts stably, so `dot` with this sort makes the same choices:
 *
 *     macOS:  clang -O2 -dynamiclib -o libstable.dylib tools/stable_qsort.c
 *             and run `dot` with DYLD_INSERT_LIBRARIES set to that file
 *     Linux:  gcc -O2 -shared -fPIC -DPRELOAD -o libstable.so tools/stable_qsort.c
 *             and run `dot` with LD_PRELOAD set to that file
 *
 * `pixi run compare --table` then finds ratios of 1.000 for edge length and area on the whole corpus.
 */
#include <stdlib.h>
#include <string.h>

#ifdef PRELOAD
void qsort(void *base, size_t n, size_t s, int (*cmp)(const void *, const void *))
#else
static void stable_qsort(void *base, size_t n, size_t s, int (*cmp)(const void *, const void *))
#endif
{
  if (n < 2)
    return;
  char *a = base, *t = malloc(n * s);
  for (size_t w = 1; w < n; w *= 2) {
    for (size_t lo = 0; lo < n; lo += 2 * w) {
      size_t mid = lo + w < n ? lo + w : n, hi = lo + 2 * w < n ? lo + 2 * w : n;
      size_t i = lo, j = mid;
      for (size_t k = lo; k < hi; k++) {
        if (i < mid && (j >= hi || cmp(a + i * s, a + j * s) <= 0)) {
          memcpy(t + k * s, a + i * s, s);
          i++;
        } else {
          memcpy(t + k * s, a + j * s, s);
          j++;
        }
      }
    }
    memcpy(a, t, n * s);
  }
  free(t);
}

#ifndef PRELOAD
__attribute__((used)) static struct {
  const void *replacement, *replaced;
} interposers[] __attribute__((section("__DATA,__interpose"))) = {
    {(const void *)stable_qsort, (const void *)qsort}};
#endif
