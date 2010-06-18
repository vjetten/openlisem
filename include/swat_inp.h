
const PROFILE *ProfileNr(int profileNr);
/* RET profile or NULL if not found  */

int NrZoneNodes(void);
/* returns nr of Nodes 
 * or -1 if soil tables  are not loaded yet
 * function only valid in 
 * the current implementation
 * where all pixels have an equal 
 * zone definition
 */
