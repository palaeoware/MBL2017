# MBL2017 - Inno Setup
## Phylogenetic tree and matrix simulation based on birth–death systems

Coded by:
 - Mark Sutton (m.sutton@imperial.ac.uk)

Additional coding:
 - GUI theme designed by Alan R.T. Spencer.
 - Tinkering, GUI implementation and release by Russell J. Garwood.

Concepts and algorithms:
- Mark Sutton
- Julia Sigwart (j.sigwart@qub.ac.uk)

______

## Relevant references:
Sigwart, J.D., Sutton, M.D. and Bennett, K.D., 2017. How big is a genus? Towards a nomothetic systematics. Zoological Journal of the Linnean Society.

Keating, J.N., Sansom, R.S., Sutton, M.D., Knight, C.G. and Garwood, R.J. 2019. Comparing methods of morphological phylogenetic estimation using novel evolutionary simulations. 
_____


CONTENTS:

1. Creating Windows Installer Program - Inno Setup

_____

## 1. Creating Windows Installer Program - Inno Setup

Inno Setup is a free installer for Windows programs by Jordan Russell and Martijn Laan. First introduced in 1997, Inno Setup today rivals and even surpasses many commercial installers in feature set and stability.

The program can be downloaded from: [http://jrsoftware.org/isinfo.php](http://jrsoftware.org/isinfo.php)

The file MBL2017installer.iss contained in this folder should be modified and used to create the Windows Installation Binary. It expects all support .DLLs and the .EXE to be placed in a ./bin folder. The Inno Setup will create the new install binary in a folder called ./build

The MBL2017installer.iss should be altered with the latest version number etc... before being run.

_____


t:@palaeoware

e:palaeoware@gmail.com

w:https://github.com/palaeoware
