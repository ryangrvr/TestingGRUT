// Exact signed permutation canonicalization; no floating-point arithmetic.
#include <cstdint>
#include <vector>

extern "C" int canonical_state(
    uint64_t bits, const int32_t* field, int width, int legs,
    const int32_t* ids, int num_ids, const int32_t* modes,
    const int32_t* legmaps, const uint64_t* flipmasks,
    const int32_t* backgrounds, int group_size,
    uint64_t* bestbits, int32_t* bestfield, int32_t* sign, int32_t* orbit) {
    if (width > 63 || num_ids <= 0 || legs <= 0) return -1;
    bool have=false, conflict=false;
    int count=0, bestsign=0;
    std::vector<int32_t> out(legs);
    for (int a=0; a<num_ids; ++a) {
        int p=ids[a];
        uint64_t transformed=0, remaining=bits;
        int parity=__builtin_popcountll(bits & flipmasks[p]);
        while (remaining) {
            uint64_t bit=remaining & -remaining;
            int j=__builtin_ctzll(bit), target=modes[p*width+j];
            parity+=__builtin_popcountll(transformed >> target);
            transformed |= (uint64_t(1)<<target);
            remaining ^= bit;
        }
        for (int j=0; j<legs; ++j) out[legmaps[p*legs+j]]=field[j];
        int phase=backgrounds[p]*((parity & 1) ? -1 : 1);
        int comparison=0;
        if (!have || transformed < *bestbits) comparison=-1;
        else if (transformed > *bestbits) comparison=1;
        else {
            for (int j=0; j<legs; ++j) {
                if (out[j]<bestfield[j]) { comparison=-1; break; }
                if (out[j]>bestfield[j]) { comparison=1; break; }
            }
        }
        if (comparison<0) {
            have=true; *bestbits=transformed;
            for (int j=0; j<legs; ++j) bestfield[j]=out[j];
            bestsign=phase; count=1; conflict=false;
        } else if (comparison==0) {
            ++count; if (phase!=bestsign) conflict=true;
        }
    }
    if (!count || group_size % count) return -2;
    *sign=conflict ? 0 : bestsign;
    *orbit=group_size/count;
    return 0;
}
