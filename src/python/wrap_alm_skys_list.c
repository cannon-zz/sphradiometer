SWIGINTERN PyObject *_wrap_generate_alm_skys(PyObject *SWIGUNUSEDPARM(self), PyObject *args) {
  PyObject *resultobj = 0;
  PyObject *result1 = 0;
  PyObject *result2 = 0;
  struct sh_series *arg1 = NULL ;
  struct sh_series *arg2 = NULL ;
  struct correlator_network_plan_fd *arg3 = (struct correlator_network_plan_fd *) 0 ;
  struct correlator_network_plan_fd *arg4 = (struct correlator_network_plan_fd *) 0 ;
  COMPLEX16TimeSeries **arg5 = (COMPLEX16TimeSeries **) 0 ;
  COMPLEX16Sequence **arg6 = (COMPLEX16Sequence **) 0 ;
  struct sh_series *arg7 = (struct sh_series *) 0 ;
  void *argp3 = 0 ;
  int res3 = 0 ;
  void *argp4 = 0 ;
  int res4 = 0 ;
  void *argp5 = 0 ;
  int res5 = 0 ;
  void *argp6 = 0 ;
  int res6 = 0 ;
  void *argp7 = 0 ;
  int res7 = 0 ;
  PyObject * obj2 = 0 ;
  PyObject * obj3 = 0 ;
  PyObject * obj4 = 0 ;
  PyObject * obj5 = 0 ;
  PyObject * obj6 = 0 ;
  int result;
  int i;

  if (!PyArg_ParseTuple(args,(char *)"OOOOO:generate_alm_skys",&obj2,&obj3,&obj4,&obj5,&obj6)) SWIG_fail;

  res3 = SWIG_ConvertPtr(obj2, &argp3,SWIGTYPE_p_correlator_network_plan_fd, 0 |  0 );
  if (!SWIG_IsOK(res3)) {
    SWIG_exception_fail(SWIG_ArgError(res3), "in method '" "generate_alm_skys" "', argument " "3"" of type '" "struct correlator_network_plan_fd *""'"); 
  }
  arg3 = (struct correlator_network_plan_fd *)(argp3);

  res4 = SWIG_ConvertPtr(obj3, &argp4,SWIGTYPE_p_correlator_network_plan_fd, 0 |  0 );
  if (!SWIG_IsOK(res4)) {
    SWIG_exception_fail(SWIG_ArgError(res4), "in method '" "generate_alm_skys" "', argument " "4"" of type '" "struct correlator_network_plan_fd *""'"); 
  }
  arg4 = (struct correlator_network_plan_fd *)(argp4);

  arg5 = malloc(PyList_Size(obj4) * sizeof(*arg5));
  for(i = 0; i < PyList_Size(obj4); i++) {
    PyObject *item = PyList_GetItem(obj4, i);
    fprintf(stderr, "AAAA\n");
    res5 = SWIG_ConvertPtr(item, &argp5, SWIGTYPE_p_COMPLEX16TimeSeries, 0 | 0);
    fprintf(stderr, "AAAA1\n");
    if (!SWIG_IsOK(res5)) {
      SWIG_exception_fail(SWIG_ArgError(res5), "in method '" "generate_alm_skys" "', argument " "5"" of type '" "list of COMPLEX16TimeSeries *""'");
    }
    arg5[i] = argp5;
  }

  arg6 = malloc(PyList_Size(obj5) * sizeof(*arg6));
  for(i = 0; i < PyList_Size(obj5); i++) {
    PyObject *item = PyList_GetItem(obj5, i);
    res6 = SWIG_ConvertPtr(item, &argp6, SWIGTYPE_p_COMPLEX16Sequence, 0 | 0);
    if (!SWIG_IsOK(res6)) {
      SWIG_exception_fail(SWIG_ArgError(res6), "in method '" "generate_alm_skys" "', argument " "6"" of type '" "list of COMPLEX16Sequence *""'");
    }
    arg6[i] = argp6;
  }

  res7 = SWIG_ConvertPtr(obj6, &argp7,SWIGTYPE_p_sh_series, 0 |  0 );
  if (!SWIG_IsOK(res7)) {
    SWIG_exception_fail(SWIG_ArgError(res7), "in method '" "generate_alm_skys" "', argument " "7"" of type '" "struct sh_series *""'"); 
  }
  arg7 = (struct sh_series *)(argp7);

  fprintf(stderr, "DDDD\n");
  result = (int)generate_alm_skys(&arg1,&arg2,arg3,arg4,arg5,arg6,arg7);
  free(arg5);
  free(arg6);
  fprintf(stderr, "EEEE\n");
  if(result) {
    /* FIXME: raise exception */
    fprintf(stderr, "error generate_alm_skys\n");
  }
  result1 = PyCapsule_New((void *)arg1, NULL, NULL);
  result2 = PyCapsule_New((void *)arg2, NULL, NULL);
  resultobj = PyTuple_Pack(2, result1, result2);
  return resultobj;
fail:
  return NULL;
}
