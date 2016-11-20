#!/bin/sh
#append $NWELLS if given
cmd=/$NWELLS/$1
shift
$cmd $@
