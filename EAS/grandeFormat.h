UShort_t linkID, packSize, channel, ndata; //  16 bit vars
UInt_t evtID, evtIDlast = -1, errorsTot, clNumber, VMEadress, end;   // 32 bit vars
ULong_t time;//64 bit vars
UShort_t * data = NULL;
int nPacksTot = 0;
int data_time[4];
UInt_t h, m, s, mls, mks, ns, Delay;
