#include "node.h"
#include <string>

int ** initializeHashMtx();

int hashVariants (string trio, char shift);

pair<string, char> unhashVariants (int variantHash);

void checkVariants(LeafData* pCurrentData);


