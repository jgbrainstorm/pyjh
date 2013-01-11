# This file was automatically generated by SWIG (http://www.swig.org).
# Version 1.3.40
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.
# This file is compatible with both classic and new-style classes.

from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_ecgmm', [dirname(__file__)])
        except ImportError:
            import _ecgmm
            return _ecgmm
        if fp is not None:
            try:
                _mod = imp.load_module('_ecgmm', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _ecgmm = swig_import_helper()
    del swig_import_helper
else:
    import _ecgmm
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


class SwigPyIterator(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, SwigPyIterator, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, SwigPyIterator, name)
    def __init__(self, *args, **kwargs): raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _ecgmm.delete_SwigPyIterator
    __del__ = lambda self : None;
    def value(self): return _ecgmm.SwigPyIterator_value(self)
    def incr(self, n = 1): return _ecgmm.SwigPyIterator_incr(self, n)
    def decr(self, n = 1): return _ecgmm.SwigPyIterator_decr(self, n)
    def distance(self, *args): return _ecgmm.SwigPyIterator_distance(self, *args)
    def equal(self, *args): return _ecgmm.SwigPyIterator_equal(self, *args)
    def copy(self): return _ecgmm.SwigPyIterator_copy(self)
    def next(self): return _ecgmm.SwigPyIterator_next(self)
    def __next__(self): return _ecgmm.SwigPyIterator___next__(self)
    def previous(self): return _ecgmm.SwigPyIterator_previous(self)
    def advance(self, *args): return _ecgmm.SwigPyIterator_advance(self, *args)
    def __eq__(self, *args): return _ecgmm.SwigPyIterator___eq__(self, *args)
    def __ne__(self, *args): return _ecgmm.SwigPyIterator___ne__(self, *args)
    def __iadd__(self, *args): return _ecgmm.SwigPyIterator___iadd__(self, *args)
    def __isub__(self, *args): return _ecgmm.SwigPyIterator___isub__(self, *args)
    def __add__(self, *args): return _ecgmm.SwigPyIterator___add__(self, *args)
    def __sub__(self, *args): return _ecgmm.SwigPyIterator___sub__(self, *args)
    def __iter__(self): return self
SwigPyIterator_swigregister = _ecgmm.SwigPyIterator_swigregister
SwigPyIterator_swigregister(SwigPyIterator)

class IntVector(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, IntVector, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, IntVector, name)
    __repr__ = _swig_repr
    def iterator(self): return _ecgmm.IntVector_iterator(self)
    def __iter__(self): return self.iterator()
    def __nonzero__(self): return _ecgmm.IntVector___nonzero__(self)
    def __bool__(self): return _ecgmm.IntVector___bool__(self)
    def __len__(self): return _ecgmm.IntVector___len__(self)
    def pop(self): return _ecgmm.IntVector_pop(self)
    def __getslice__(self, *args): return _ecgmm.IntVector___getslice__(self, *args)
    def __setslice__(self, *args): return _ecgmm.IntVector___setslice__(self, *args)
    def __delslice__(self, *args): return _ecgmm.IntVector___delslice__(self, *args)
    def __delitem__(self, *args): return _ecgmm.IntVector___delitem__(self, *args)
    def __getitem__(self, *args): return _ecgmm.IntVector___getitem__(self, *args)
    def __setitem__(self, *args): return _ecgmm.IntVector___setitem__(self, *args)
    def append(self, *args): return _ecgmm.IntVector_append(self, *args)
    def empty(self): return _ecgmm.IntVector_empty(self)
    def size(self): return _ecgmm.IntVector_size(self)
    def clear(self): return _ecgmm.IntVector_clear(self)
    def swap(self, *args): return _ecgmm.IntVector_swap(self, *args)
    def get_allocator(self): return _ecgmm.IntVector_get_allocator(self)
    def begin(self): return _ecgmm.IntVector_begin(self)
    def end(self): return _ecgmm.IntVector_end(self)
    def rbegin(self): return _ecgmm.IntVector_rbegin(self)
    def rend(self): return _ecgmm.IntVector_rend(self)
    def pop_back(self): return _ecgmm.IntVector_pop_back(self)
    def erase(self, *args): return _ecgmm.IntVector_erase(self, *args)
    def __init__(self, *args): 
        this = _ecgmm.new_IntVector(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(self, *args): return _ecgmm.IntVector_push_back(self, *args)
    def front(self): return _ecgmm.IntVector_front(self)
    def back(self): return _ecgmm.IntVector_back(self)
    def assign(self, *args): return _ecgmm.IntVector_assign(self, *args)
    def resize(self, *args): return _ecgmm.IntVector_resize(self, *args)
    def insert(self, *args): return _ecgmm.IntVector_insert(self, *args)
    def reserve(self, *args): return _ecgmm.IntVector_reserve(self, *args)
    def capacity(self): return _ecgmm.IntVector_capacity(self)
    __swig_destroy__ = _ecgmm.delete_IntVector
    __del__ = lambda self : None;
IntVector_swigregister = _ecgmm.IntVector_swigregister
IntVector_swigregister(IntVector)

class DoubleVector(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, DoubleVector, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, DoubleVector, name)
    __repr__ = _swig_repr
    def iterator(self): return _ecgmm.DoubleVector_iterator(self)
    def __iter__(self): return self.iterator()
    def __nonzero__(self): return _ecgmm.DoubleVector___nonzero__(self)
    def __bool__(self): return _ecgmm.DoubleVector___bool__(self)
    def __len__(self): return _ecgmm.DoubleVector___len__(self)
    def pop(self): return _ecgmm.DoubleVector_pop(self)
    def __getslice__(self, *args): return _ecgmm.DoubleVector___getslice__(self, *args)
    def __setslice__(self, *args): return _ecgmm.DoubleVector___setslice__(self, *args)
    def __delslice__(self, *args): return _ecgmm.DoubleVector___delslice__(self, *args)
    def __delitem__(self, *args): return _ecgmm.DoubleVector___delitem__(self, *args)
    def __getitem__(self, *args): return _ecgmm.DoubleVector___getitem__(self, *args)
    def __setitem__(self, *args): return _ecgmm.DoubleVector___setitem__(self, *args)
    def append(self, *args): return _ecgmm.DoubleVector_append(self, *args)
    def empty(self): return _ecgmm.DoubleVector_empty(self)
    def size(self): return _ecgmm.DoubleVector_size(self)
    def clear(self): return _ecgmm.DoubleVector_clear(self)
    def swap(self, *args): return _ecgmm.DoubleVector_swap(self, *args)
    def get_allocator(self): return _ecgmm.DoubleVector_get_allocator(self)
    def begin(self): return _ecgmm.DoubleVector_begin(self)
    def end(self): return _ecgmm.DoubleVector_end(self)
    def rbegin(self): return _ecgmm.DoubleVector_rbegin(self)
    def rend(self): return _ecgmm.DoubleVector_rend(self)
    def pop_back(self): return _ecgmm.DoubleVector_pop_back(self)
    def erase(self, *args): return _ecgmm.DoubleVector_erase(self, *args)
    def __init__(self, *args): 
        this = _ecgmm.new_DoubleVector(*args)
        try: self.this.append(this)
        except: self.this = this
    def push_back(self, *args): return _ecgmm.DoubleVector_push_back(self, *args)
    def front(self): return _ecgmm.DoubleVector_front(self)
    def back(self): return _ecgmm.DoubleVector_back(self)
    def assign(self, *args): return _ecgmm.DoubleVector_assign(self, *args)
    def resize(self, *args): return _ecgmm.DoubleVector_resize(self, *args)
    def insert(self, *args): return _ecgmm.DoubleVector_insert(self, *args)
    def reserve(self, *args): return _ecgmm.DoubleVector_reserve(self, *args)
    def capacity(self): return _ecgmm.DoubleVector_capacity(self)
    __swig_destroy__ = _ecgmm.delete_DoubleVector
    __del__ = lambda self : None;
DoubleVector_swigregister = _ecgmm.DoubleVector_swigregister
DoubleVector_swigregister(DoubleVector)


def BICecgmm(*args):
  return _ecgmm.BICecgmm(*args)
BICecgmm = _ecgmm.BICecgmm

def AICecgmm(*args):
  return _ecgmm.AICecgmm(*args)
AICecgmm = _ecgmm.AICecgmm

def BICecgmmFixedComponent(*args):
  return _ecgmm.BICecgmmFixedComponent(*args)
BICecgmmFixedComponent = _ecgmm.BICecgmmFixedComponent

def AICecgmmFixedComponent(*args):
  return _ecgmm.AICecgmmFixedComponent(*args)
AICecgmmFixedComponent = _ecgmm.AICecgmmFixedComponent


