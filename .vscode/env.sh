#!/bin/bash

# Load .bashrc like a normal shell

[ -f ~/.bashrc ] && . ~/.bashrc

# Source .vscode/env.local.sh if it exists

_local="$(dirname "${BASH_SOURCE[0]}")/env.local.sh"

[ -f "$_local" ] && . "$_local"

unset _local
