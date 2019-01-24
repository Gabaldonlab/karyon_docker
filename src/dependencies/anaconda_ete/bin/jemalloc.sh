#!/bin/sh

prefix=/nfs/users/tg/mnaranjo/anaconda_ete
exec_prefix=/nfs/users/tg/mnaranjo/anaconda_ete
libdir=${exec_prefix}/lib

LD_PRELOAD=${libdir}/libjemalloc.so.2
export LD_PRELOAD
exec "$@"
