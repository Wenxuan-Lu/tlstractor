AC_INIT([tlstractor], [0.1.0], [wenxuanlu@hotmail.com])
AC_CONFIG_SRCDIR([DESCRIPTION])
AC_CONFIG_FILES([src/Makevars])
AC_PROG_CC
AC_LANG([C])

: ${ZLIB_CPPFLAGS=""}
: ${ZLIB_LIBS=""}
TLSTRACTOR_CPPFLAGS=""
TLSTRACTOR_LIBS=""

dnl Allow manual override, e.g.
dnl R CMD INSTALL . --configure-vars='ZLIB_CPPFLAGS=-I/opt/include ZLIB_LIBS=-L/opt/lib -lz'
AC_ARG_VAR([ZLIB_CPPFLAGS], [C preprocessor flags for zlib headers (for example -I/path/include)])
AC_ARG_VAR([ZLIB_LIBS], [linker flags for zlib (for example -L/path/lib -lz)])

if test -n "$ZLIB_CPPFLAGS" || test -n "$ZLIB_LIBS"; then
  AC_MSG_NOTICE([using user-provided zlib flags from --configure-vars or environment (any missing part may still be filled by pkg-config)])
fi

dnl Prefer pkg-config when available
AC_PATH_PROG([PKG_CONFIG], [pkg-config], [], [$PATH:/usr/bin:/usr/local/bin:/opt/bin:/opt/homebrew/bin])

dnl If either include flags or libs are missing, try pkg-config to fill only the missing pieces.
if test -z "$ZLIB_CPPFLAGS" || test -z "$ZLIB_LIBS"; then
    if test -n "$PKG_CONFIG"; then
        AC_MSG_CHECKING([whether pkg-config can find zlib])
        if "$PKG_CONFIG" --exists zlib; then
            AC_MSG_RESULT([yes])
            if test -z "$ZLIB_CPPFLAGS"; then
                ZLIB_CPPFLAGS=`"$PKG_CONFIG" --cflags zlib`
            fi
            if test -z "$ZLIB_LIBS"; then
                ZLIB_LIBS=`"$PKG_CONFIG" --libs zlib`
            fi
        else
            AC_MSG_RESULT([no])
        fi
    fi
fi

CPPFLAGS_SAVED="$CPPFLAGS"
LIBS_SAVED="$LIBS"

CPPFLAGS="$CPPFLAGS $ZLIB_CPPFLAGS"
LIBS="$LIBS $ZLIB_LIBS"

AC_CHECK_HEADER([zlib.h], [],
    [AC_MSG_ERROR([
zlib headers not found (missing zlib.h).

What this means:
  - zlib development headers are not installed, or
  - the include path for zlib is not visible to the compiler.

How to resolve:
  - macOS (Homebrew): brew install zlib pkg-config
    then reinstall with:
    R CMD INSTALL . --configure-vars='ZLIB_CPPFLAGS=-I/opt/homebrew/opt/zlib/include ZLIB_LIBS=-L/opt/homebrew/opt/zlib/lib -lz'
  - Debian/Ubuntu: sudo apt-get install zlib1g-dev pkg-config
  - Fedora/RHEL/CentOS: sudo dnf install zlib-devel pkgconf-pkg-config
  - Alpine: apk add zlib-dev pkgconf
  - Conda: conda install -c conda-forge zlib pkg-config

If zlib is already installed in a custom location, provide include/library flags with --configure-vars.
Examples:
  R CMD INSTALL . --configure-vars='ZLIB_CPPFLAGS=-I/path/include ZLIB_LIBS=-L/path/lib -lz'
  R CMD INSTALL . --configure-vars='ZLIB_LIBS=-L/path/lib -lz'   (leave CPPFLAGS to auto-detect)
  ZLIB_CPPFLAGS='-I/path/include' ZLIB_LIBS='-L/path/lib -lz' R CMD INSTALL .
])])

AC_SEARCH_LIBS([deflate], [z], [],
    [AC_MSG_ERROR([
zlib library not found (failed to link symbol deflate, typically from -lz).

What this means:
  - the zlib development library is not installed, or
  - the library path for zlib is not visible to the linker.

How to resolve:
  - macOS (Homebrew): brew install zlib pkg-config
  - Debian/Ubuntu: sudo apt-get install zlib1g-dev pkg-config
  - Fedora/RHEL/CentOS: sudo dnf install zlib-devel pkgconf-pkg-config
  - Alpine: apk add zlib-dev pkgconf
  - Conda: conda install -c conda-forge zlib pkg-config

Then reinstall with manual flags if needed:
  R CMD INSTALL . --configure-vars='ZLIB_CPPFLAGS=-I/path/include ZLIB_LIBS=-L/path/lib -lz'
  R CMD INSTALL . --configure-vars='ZLIB_LIBS=-L/path/lib -lz'   (if headers are already found)
  ZLIB_CPPFLAGS='-I/path/include' ZLIB_LIBS='-L/path/lib -lz' R CMD INSTALL .

Note for Windows:
  - This configure script is used on Unix-like platforms.
  - Windows package installation may require configure.win and/or src/Makevars.win.
])])

if test -z "$ZLIB_LIBS"; then
    case "$ac_cv_search_deflate" in
        no) ;;
        none\ required) ZLIB_LIBS="" ;;
        *) ZLIB_LIBS="$ac_cv_search_deflate" ;;
    esac
fi

CPPFLAGS="$CPPFLAGS_SAVED"
LIBS="$LIBS_SAVED"

TLSTRACTOR_CPPFLAGS="$ZLIB_CPPFLAGS"
TLSTRACTOR_LIBS="$ZLIB_LIBS"

AC_SUBST([TLSTRACTOR_CPPFLAGS])
AC_SUBST([TLSTRACTOR_LIBS])

AC_OUTPUT
