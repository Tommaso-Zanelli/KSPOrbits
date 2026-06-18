PYTHON_INCLUDE = C:\Users\t.zanelli\AppData\Local\Programs\Python\Python313\include
PYTHON_LIBS    = C:\Users\t.zanelli\AppData\Local\Programs\Python\Python313\libs
PYBIND_INCLUDE = C:\Users\t.zanelli\AppData\Local\Programs\Python\Python313\Lib\site-packages\pybind11\include

CFLAGS = /O2 /LD /EHsc /MD
INCLUDES = /I "$(PYTHON_INCLUDE)" /I "$(PYBIND_INCLUDE)"
LDFLAGS = /link /LIBPATH:"$(PYTHON_LIBS)"

hello.pyd: hello.cpp
	cl $(CFLAGS) $(INCLUDES) hello.cpp $(LDFLAGS) /OUT:hello.pyd

clean:
    del hello.exp hello.lib hello.obj hello.pyd
