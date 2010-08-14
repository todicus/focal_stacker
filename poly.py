# Implementation of polynomials

class Error(Exception):
    pass

class Polynomial:
    "An immutable polynomial class"
    def __init__(self, seq):
        self._coeffs = tuple(seq)

    def __len__(self):
        return len(self._coeffs)

    def getcoeff(self, key):
        try:
            return self._coeffs[key]
        except IndexError:
            return 0.0

    def __iter__(self):
        return iter(self._coeffs)

    def __getitem__(self, key):
        if isinstance(key, slice):
            if key.start < 0 or key.stop < 0:
                raise Error('Negative indicies not supported by Polynomial')
            start = key.start
            stop = key.stop
            if stop == sys.maxint:
                stop = len(self._coeffs)
            step = key.step
            if step is None:
                step = 1
            if step < 0:
                raise Error('Negative steps not supported by Polynomial')
            indices = range(start, stop, step)
            return [self.getcoeff(x) for x in indices]
        else:
            return self.getcoeff(key)

    def perform2(self, op, other):
        "Perform an element-wise operation on two operands (self, other)"
        outlen = max(len(self),len(other))
        return Polynomial([op(self[i],other[i]) for i in range(outlen)])

    def perform1(self, op):
        "Perform an element-wise operation on one operand (self)"
        return Polynomial([op(self[i]) for i in range(len(self))])

    def __add__(self, other):
        "Add Polynomials"
        return self.perform2(lambda x,y: x+y, other)

    def __sub__(self, other):
        "Subtract polynomials"
        return self.perform2(lambda x,y: x-y, other)

    def __radd__(self, other):
        "Add Polynomials"
        return self.perform2(lambda x,y: y+x, other)

    def __rsub__(self, other):
        "Subtract polynomials"
        return self.perform2(lambda x,y: y-x, other)

    def __mul__(self, other):
        "Scale polynomial"
        return self.perform1(lambda x: other*x)

    def __rmul__(self, other):
        "Scale polynomial"
        return self.perform1(lambda x: other*x)

    def __lshift__(self, other):
        "Shift polynomial to a higher degree"
        return Polynomial((0,)*other + self._coeffs)

    def __rshift__(self, other):
        "Shift polynomial to a lower degree (truncate)"
        return Polynomial(self[other:])

    def __str__(self):
        "Display the Polynomial"
        return 'Polynomial['+','.join([('%g' % self[i]) for i in range(len(self))])+']'

    def __repr__(self):
        "Display the Polynomial"
        return 'Polynomial['+repr(self._coeffs)+']'

