# MLKMC - Machine-learning augmented kinetic Monte Carlo
# 
# Copyright (C) 2026 Jyri Kimari
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

include MAKE/makefile.defs

all: ubuntu

ubuntu: 
	make -f MAKE/makefile.ubuntu

puhti: 
	make -f MAKE/makefile.puhti

tetralith:
	make -f MAKE/makefile.tetralith

turso:
	make -f MAKE/makefile.turso

clean:
	rm *.o *.mod

clean-all:
	rm ${TARGETS} *.o *.mod
