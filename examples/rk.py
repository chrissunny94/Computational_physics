class RK4(object):

    def __init__(self, *functions):
        self.f = functions
        self.t = 0

    def solve(self, y, h, n):
        t = []
        res = [[] for _ in y]

        while self.t <= n and h != 0:
            t.append(self.t)
            y = self._solve(y, self.t, h)

            for c, val in enumerate(y):
                res[c].append(val)

            self.t += h
            if self.t + h > n:
                h = n - self.t

        return t, res

    def _solve(self, y, t, h):
        f = self.f

        # k1
        k1 = [h * fn(t, *y) for fn in f]

        # k2
        k2 = [
            h * fn(
                t + 0.5 * h,
                *[y[i] + 0.5 * k1[i] for i in range(len(y))]
            )
            for fn in f
        ]

        # k3
        k3 = [
            h * fn(
                t + 0.5 * h,
                *[y[i] + 0.5 * k2[i] for i in range(len(y))]
            )
            for fn in f
        ]

        # k4
        k4 = [
            h * fn(
                t + h,
                *[y[i] + k3[i] for i in range(len(y))]
            )
            for fn in f
        ]

        # update
        return [
            y[i] + (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) / 6.0
            for i in range(len(y))
        ]
