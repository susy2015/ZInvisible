AC_INIT(program, 1.0)

dnl Switch to a C++ compiler, and check if it works.
AC_LANG(C++)
AC_PROG_CXX

dnl set --std=c++11
case "$CXXFLAGS" in 
  *"--std=c++11"*)
    # Do nothing
    ;;
  *)
    CXXFLAGS="$CXXFLAGS --std=c++11" 
esac

dnl check for ldconfig
AC_CHECK_PROG(ldconfigfound, ldconfig, yes, no)

LDCONFIG=
if test ${ldconfigfound} == yes; then
   LDCONFIG="ldconfig"
fi

AC_SUBST(LDCONFIG, $LDCONFIG)

dnl check if this is a mac
LIBSUFFIX="so"
OSNAME=$(uname -s)
if [[ "$OSNAME" == "Darwin" ]]; then
   LIBSUFFIX="dylib"
fi

AC_SUBST(LIBSUFFIX, $LIBSUFFIX)

dnl check for root
AC_CHECK_PROG(rootconfigfound, root-config, yes, no)
AC_CHECK_PROG(rootcintfound, rootcint, yes, no)

if test ${rootconfigfound} == no; then
   AC_ERROR([ROOT: root-config not found.])
fi

if test ${rootcintfound} == no; then
   AC_ERROR([ROOT: rootcint not found.])
fi   

dnl check if scram is avaliable, and if so, get the numpy include path
AC_CHECK_PROG(scramfound, scram, yes, no)

if test ${scramfound} == yes; then
   AC_SUBST(SCRAMNUMPYCFG, $(scram tool tag py2-numpy PY2_NUMPY_BASE))
else
   AC_SUBST(SCRAMNUMPYCFG, "")
fi

dnl check for python 
AC_CHECK_PROG(pythonconfigfound, python-config, yes, no)

if test ${pythonconfigfound} == no; then
dnl   AC_WARN([Optional dependency missing: python-config not found.])
   AC_SUBST(PYTHONCFG, "")
else
   AC_SUBST(PYTHONCFG, python-config)
fi

dnl check for top tagger code
if test -z ${TOPTAGGERDIR+x}; then
   TOPTAGGERDIR="${PWD}/../../TopTagger"
fi

if test -z ${ZINVISIBLE}; then 
    ZINVISIBLE="${PWD}"
fi

AC_CHECK_HEADER($TOPTAGGERDIR/TopTagger/interface/TopTagger.h,        [], [AC_ERROR([TopTagger includes are not found!!!!!!!!!!])])
AC_CHECK_HEADERS($TOPTAGGERDIR/CfgParser/include/TTException.h,     [], [AC_ERROR([CfgParser includes are not found!!!!!!!!!!])])

AC_SUBST(TAGGERDIR, $TOPTAGGERDIR)

dnl check for SusyAnaTools
if test -z ${SUSYANATOOLSDIR+x}; then
   SUSYANATOOLSDIR="${PWD}/../../SusyAnaTools"
fi

AC_CHECK_HEADER($SUSYANATOOLSDIR/Tools/samples.h,        [], [AC_ERROR([SusyAnaTools includes are not found!!!!!!!!!!])])
AC_CHECK_FILE($SUSYANATOOLSDIR/Tools/NTupleReader.h,   [], [AC_ERROR([SusyAnaTools includes are not found!!!!!!!!!!])])

AC_SUBST(SUSYANATOOLSDIR, $SUSYANATOOLSDIR)

dnl check for LHAPDF

AC_CHECK_FILE($LHAPDF_DATA_PATH/../../include/LHAPDF/LHAPDF.h, [], [AC_ERROR([LHAPDF includes are not found!!!!!!!!!!])])

AC_SUBST(LHAPDFDIR, $LHAPDF_DATA_PATH)

dnl Process Makefile.in to create Makefile
AC_OUTPUT(Makefile)
