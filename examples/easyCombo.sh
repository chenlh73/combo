#!/bin/bash
# ============================================================
# easyCombo.sh — qcombo commutator calculation script
#
# Usage: bash easyCombo.sh
# ============================================================

# qcombo command format:
#   qcombo <left_body> <right_body> [options]
#
# Options:
#   -c, --contraction  Contraction body: single int(0), range(0-2), comma-sep(0,1,2), all
#   -w, --wick-mode    Wick mode: MR(multi-ref, default) / SR(single-ref)
#   -p, --parallel          Enable parallel computing (recommended for body rank >= 3)
#   -ns, --no-process       Disable progress bar
#   -q, --quiet        Quiet mode
#   -lo, --latex-output     Custom .tex output path
#   -ao, --amc-output       Custom .amc output path
#   -i, --interactive  Interactive mode
#   -v, --version      Version info
#   -h, --help         Help info

# ===================== Example tasks =====================

# [1,1] — 1-body × 1-body, simplest usage
qcombo 1 1

# Full parameter input
# qcombo 1 1 -c 0-1 -lo commutator_11.tex -ao commutator_11.amc

# More display options
# qcombo 1 1 -c 0,1 -lo commutator_11.tex -ao commutator_11.amc -w MR

# Enable parallel computing for large body rank calculations
# qcombo 2 3 --p

# Interactive mode
# qcombo -i
# or
# qcombo
