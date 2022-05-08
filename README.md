# jstar-lbvs.cpp
[jstar]'s daemon for ligand-based virtual screening (LBVS), written in C++.

A JavaScript implementation is available at https://github.com/HongjianLi/jstar-lbvs.js

## Usage
    bin/lbvs ../jstar/databases localhost 27017 jstard password

## Architecture
![jstar architecture](https://github.com/HongjianLi/jstar/blob/master/public/architecture.png)

## Components
### Database
* [MongoDB]
### Daemon
* [MongoDB C++ Driver]
* [rdkit]
* [boost]

## License
[MIT License]

## Developers
* [Jacky Lee]
* Jessie Sze

## Logo
![jstar logo](https://github.com/HongjianLi/jstar/blob/master/public/logo.svg)

[jstar]: https://github.com/HongjianLi/jstar
[MongoDB]: https://github.com/mongodb/mongo
[MongoDB C++ Driver]: https://github.com/mongodb/mongo-cxx-driver
[rdkit]: https://github.com/rdkit/rdkit
[boost]: https://github.com/boostorg/boost
[MIT License]: https://github.com/HongjianLi/jstar-lbvs.cpp/blob/master/LICENSE
[Jacky Lee]: https://github.com/HongjianLi
