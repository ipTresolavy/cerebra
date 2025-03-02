#!/bin/bash

SCRIPTPATH="$(
  cd -- "$(dirname "$0")" >/dev/null 2>&1 || return 255
  pwd -P
)"

find $SCRIPTPATH -not \( -path "*build*" -prune \) -not \( -path "*bin*" -prune \) -name "*.c" -o -name "*.h" -type f | xargs clang-format -i --style=file --verbose
