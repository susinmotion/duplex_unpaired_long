#include "node.h"
#include "variant.h"
#include <string>
#include <vector>
int ** initializeHashMtx();

int hashVariants (string trio, char shift);

pair<string, char> unhashVariants (int variantHash);

void checkVariants(LeafData* pCurrentData);

vector  <Variant*> bowtieCheckVariants(string sequence, string gene="p53");
