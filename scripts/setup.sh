#!/bin/sh
#append $NWELLS if given
mv $NWELLS/* /usr/local/bin
#run the command
"$@"
