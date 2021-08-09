#ifndef BASE_H
#define BASE_H

namespace SS {
    enum IDLevel {
        IDdefault = -1,
        IDveto = 0, // for Z-veto
        IDfakableNoIso = 1,
        IDfakable = 2, // for fake background + jet cleaning
        IDtightNoIso = 3,
        IDtight = 4 // for analysis
    };
}

#endif
