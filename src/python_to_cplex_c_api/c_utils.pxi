from libc.stdlib cimport malloc, free

cdef class ArrayOfDouble():
    cdef double* impl
    cdef Py_ssize_t size
    cdef bint owner

    def __cinit__(self, values = None):

        if values is None:
            self.owner = False
            self.impl = NULL
            self.size = 0
            self.owner = False
            return

        self.owner = True
        
        if isinstance(values, int):
            self.size = values
            self.impl = <double*> malloc(self.size * sizeof(double))
            return

        self.size = len(values) if values is not None else 0
        
        if self.size == 0: 
            self.impl = NULL
            return
        
        self.impl = <double*> malloc(self.size * sizeof(double))
        for i in range(self.size): self.impl[i] = values[i]
    
    def __dealloc__(self):
        if self.impl != NULL and self.owner is True:
            free(self.impl)

    def set_size(self, size):
        self.size = size
    
    def to_list(self):
        return [self.impl[i] for i in range(self.size)]
    
    @staticmethod
    cdef from_ptr(double* impl, Py_ssize_t size):
        cdef ArrayOfDouble result = ArrayOfDouble.__new__(ArrayOfDouble)
        result.impl = impl
        result.size = size
        return result

cdef class ArrayOfInt():
    cdef int* impl
    cdef Py_ssize_t size

    def __cinit__(self, values):

        if isinstance(values, int):
            self.size = values
            self.impl = <int*> malloc(self.size * sizeof(int))
            return

        self.size = len(values) if values is not None else 0
        
        if self.size == 0: 
            self.impl = NULL
            return
        
        self.impl = <int*> malloc(self.size * sizeof(int))
        for i in range(self.size): self.impl[i] = values[i]
    
    def __dealloc__(self):
        if self.impl != NULL:
            free(self.impl)
    
    def to_list(self):
        return [self.impl[i] for i in range(self.size)]

cdef class ArrayOfChar():
    cdef char* impl
    cdef  Py_ssize_t size
    cdef bytes _xvalues

    def __cinit__(self, values):

        if isinstance(values, int):
            self.size = values
            self.impl = <char*> malloc(self.size * sizeof(char))
            return

        self.size = len(values) if values is not None else 0
        
        if self.size == 0:
            self.impl = NULL
            return

        if isinstance(values, list):
            self._xvalues = bytes("".join(values), "ascii")
        elif isinstance(values, str):
            self._xvalues = values.encode("ascii")

        self.impl = self._xvalues
    
    def to_list(self):
        return [chr(self.impl[i]) for i in range(self.size)]

cdef class ArrayOfString:
    cdef char** impl
    cdef int size
    cdef list _chars

    def __cinit__(self, values):

        self.size = len(values) if values is not None else 0
        
        if self.size == 0:
            self.impl = NULL
            return

        self.impl = <char**>malloc(self.size * sizeof(char*))

        self._chars = []

        for i in range(self.size):
            c = ArrayOfChar(values[i])
            self._chars.append(c)
            self.impl[i] = <char*>c.impl
          
    def __dealloc__(self):
        if self.impl != NULL:
            free(self.impl)

cdef class VoidPointer:
    cdef void* impl

    def __cinit__(self):
        self.impl = NULL

    @staticmethod
    cdef from_ptr(void* impl):
        cdef VoidPointer result = VoidPointer.__new__(VoidPointer)
        result.impl = impl
        return result
