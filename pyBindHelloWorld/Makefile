PYTHON_DIR     = C:\Users\t.zanelli\AppData\Local\Programs\Python\Python313
PYTHON_INCLUDE = $(PYTHON_DIR)\include
PYTHON_LIBS    = $(PYTHON_DIR)\libs
PYBIND_INCLUDE = $(PYTHON_DIR)\Lib\site-packages\pybind11\include

CFLAGS = /O2 /LD /EHsc /MD
INCLUDES = /I "$(PYTHON_INCLUDE)" /I "$(PYBIND_INCLUDE)"
LDFLAGS = /link /LIBPATH:"$(PYTHON_LIBS)"

hello.pyd: hello.cpp
	cl $(CFLAGS) $(INCLUDES) hello.cpp $(LDFLAGS) /OUT:hello.pyd

clean:
    del hello.exp hello.lib hello.obj hello.pyd
