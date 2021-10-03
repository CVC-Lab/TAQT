#ifndef CRITICALPOINT_H
#define CRITICALPOINT_H

class CriticalPoint {
public:
    CriticalPoint(int _id = 0, float _val = 0) {
        id  = _id;
        val = _val;
		LS = 1;
		US = 1;
		dbe = 0;
    }

	/**
	 * Copy constructor
	 */
	CriticalPoint(const CriticalPoint& cp) {
		id = cp.id;
		val = cp.val;
		LS = cp.LS;
		US = cp.US;
		dbe = cp.dbe;
	}

    inline static int compare(const void* p1, const void* p2);

    inline bool operator < (const CriticalPoint& pnt) const;

	float val;
    int id;
	int LS;			// Euler number of val-e
	int US;			// Euler number of val+e
	int dbe;		// Difference between the boundary edges of L(val-e) and L(val+e) 
};

int CriticalPoint::compare(const void* p1, const void* p2)
{
	CriticalPoint *pnt1, *pnt2;
	pnt1 = (CriticalPoint *)p1;
	pnt2 = (CriticalPoint *)p2;
	if(pnt1->val < pnt2->val) return -1;
	else if(pnt1->val > pnt2->val) return 1;
    else if(pnt1->id < pnt2->id) return -1;
    else if(pnt1->id > pnt2->id) return 1;
	return 0;
}

bool CriticalPoint::operator < (const CriticalPoint& pnt) const
{
    if(val < pnt.val) return true;
    else if(val > pnt.val) return false;
    else if(id < pnt.id) return true;
    return false;
}

#endif //CRITICALPOINT_H



