#ifndef DUAL_GRAPH_H
#define DUAL_GRAPH_H

#include <math.h>
#include <assert.h>

#include "cut.h"
#include "actree.h"
#include "volume3dcritical.h"
#include "moment.h"
#include "par/dynarray.h"

#include <vector>
#include <iostream>

#define _MOMENTS

using namespace std;

// global varibles
extern float sim_w[6];

/**
 *	A node for the dual graph.
 */
class DualNode {
public:
	/**
	 * Default Constructor
	 */
	DualNode(int _lv = 0) {
		level = _lv;
		fmax = all_min;
		fmin = all_max;
		volume = -1;
		ncrit = 0;
		fint = 0;
		parent = NULL;
		match = NULL;
		p_below = new vector<DualNode *>;
		p_above = new vector<DualNode *>;
		pruned = false;
	}

	~DualNode() {
		delete p_below;
		delete p_above;
	}
	 
	/**
	 * Merge a lower level dual node to a higher level node
	 * @param direct: 0 -- merge from under, 1 -- merge from above
	 */
	void merge(DualNode* p_node, int direct) {
		setFuncMin(p_node->fmin);
		setFuncMax(p_node->fmax);
		volume = ((volume > 0)? volume:0) + p_node->volume;
		fint += p_node->fint;
		ncrit += p_node->ncrit;
		p_node->parent = this;
		match = NULL;
		switch(direct) {
		case 0:
			bl += p_node->bl;
			break;
		case 1:
			bu += p_node->bu;
			break;
		default:
			assert(0);
		}
#ifdef _MOMENTS		
		moments = moments + p_node->moments;
#endif
	}
	
	/**
	 *	Add a neighbor of the lower range
	 */
	void addLowerNeighbor(DualNode* node) {
		//ngb_ls.push_back(node);
		p_below->push_back(node);
	}

	/**
	 *	Add a neighbor of the upper range.
	 */
	void addUpperNeighbor(DualNode* node) {
		//ngb_us.push_back(node);
		p_above->push_back(node);
	}

	/**
	 *	The number of neighbors in the lower range.
	 */
	int	lowerDegree() const {
		return (int) (p_below->size());
	}

	/**
	 *	The number of neighbors in the upper range.
	 */
	int upperDegree() const {
		return int (p_above->size());
	}

	/**
	 * Get the ith neighbor in the lower range.
	 */
	DualNode * lowerNeighbor(int i) const {
		return (*p_below)[i];
	}

	DualNode * lowerNeighbor(int i) {
		return (*p_below)[i];
	}

	/**
	 * Get the ith neighbor in the upper range.
	 */
	DualNode * upperNeighbor(int i) const {
		return (*p_above)[i];
	}

	DualNode * upperNeighbor(int i) {
		return (*p_above)[i];
	}

	/*
	 *	Get the importance of this node 
	 */
	float weight() {
		//return 0.5*(fmax-fmin) + 0.5*volume;
		return 0.1f*(fmax-fmin) + 0.9f*volume;
	}

	/*
	 *	
	 */
	void setFuncMax(float x) {
		if (x > fmax) {
			fmax = x;
		}
	}
	/*
	 *	
	 */
	void setFuncMin(float x) {
		if(x < fmin) {
			fmin = x;
		}
	}

	/*
	 *	Matchup this node with another node
	 */
	void setMatchNode(DualNode* dualnode) {
		match = dualnode;
		dualnode->match = this;
	}
	
	/**
	 * Check if this node has a match node
	 */
	bool isMatched() {
		return (match != NULL);
	}
	
	/**
	 * Check if the parent of this node has a match node
	 */
	bool isParentMatched(DualNode* dualnode) {
		if(parent == NULL) return true;	// top level node
		//return (parent == NULL || parent->match != NULL);
		return (parent->match != NULL && parent->match == dualnode->parent);
	}

	bool isVolumeComputed() {
		return (volume >= 0);
	}

	/**
	 * Connect a lower noder with an upper node
	 */
	static void connect(DualNode* node1, DualNode *node2) {
		node1->addUpperNeighbor(node2);
		node2->addLowerNeighbor(node1);
	}

	/**
	 *	The score of similarity of two dual nodes
	 */
	float simScore(DualNode* dualnode) {
		double s;

		double vsim, fsim, tsim;
		if(MAX(volume, dualnode->volume) == 0) {
			vsim = 0;
		} else {
			vsim = MIN(volume, dualnode->volume) / MAX(volume, dualnode->volume);
		}
		float tmax = MAX(fmax, dualnode->fmax);
		float tmin = MIN(fmin, dualnode->fmin);
		if(tmax == tmin) {
			fsim = 0.5;
		} else {
			fsim = MIN(fmax-fmin, dualnode->fmax-dualnode->fmin) / (tmax-tmin);
		}
		int ncmax = MAX(ncrit, dualnode->ncrit);
		int ncmin = MIN(ncrit, dualnode->ncrit);
		if(ncmax == 0) tsim = 1;
		else tsim = ncmin / (float) ncmax;
		tsim += (bl.simil(dualnode->bl) + bu.simil(dualnode->bu)) / 2;
		tsim /= 2;
		s = 0.8 * vsim + 0.0 * fsim + 0.2 * tsim;
		return float(s);
	}

	float simScore3(DualNode* dualnode) {
		float f1 = (volume == 0)? 0:(fint/volume);
		float f2 = (dualnode->volume == 0)? 0:(dualnode->fint/dualnode->volume);
		float len = MAX(1, MAX(fabs(f1), fabs(f2)));
		//printf("sim score2 : f1 = %f, f2 = %f\n", f1, f2);
		return simScore(dualnode)*(1-fabs(f1-f2)/len);
	}
	
	float simScore2(DualNode *dualnode) {
		float vsim, fsim, tsim;
		if(MAX(volume, dualnode->volume) == 0) {
			vsim = 0;
		} else {
			vsim = MIN(volume, dualnode->volume) / MAX(volume, dualnode->volume);
		}
		float tmax = MAX(fmax, dualnode->fmax);
		float tmin = MIN(fmin, dualnode->fmin);
		if(tmax == tmin) {
			fsim = 0.5;
		} else {
			fsim = MIN(fmax-fmin, dualnode->fmax-dualnode->fmin) / (tmax-tmin);
		}
		int ncmax = MAX(ncrit, dualnode->ncrit);
		int ncmin = MIN(ncrit, dualnode->ncrit);
		if(ncmax == 0) tsim = 1;
		else tsim = ncmin / (float) ncmax;
		tsim += (bl.simil(dualnode->bl) + bu.simil(dualnode->bu)) / 2;
		tsim /= 2;
	
		MomtAttrib attr1 = moments.toAttributes();
		MomtAttrib attr2 = dualnode->moments.toAttributes();
		float I1minusI2 = MAX(fabs(attr1.I1-attr2.I1), fabs(attr1.I2-attr2.I2));
		I1minusI2 = MAX(I1minusI2, fabs(attr1.I3-attr2.I3));
		float I1norm = (fabs(attr1.I1)+fabs(attr1.I2)+fabs(attr1.I3))/3;		
		float I2norm = (fabs(attr2.I1)+fabs(attr2.I2)+fabs(attr2.I3))/3;
		float isim = 1 - I1minusI2 / MAX(I1norm, I2norm);
		float rsim = 1 - 0.5f*(fabs(attr1.Fmin-attr2.Fmin)/ MAX(1, MAX(fabs(attr1.Fmin), fabs(attr2.Fmin)))
							+ fabs(attr1.Fmax-attr2.Fmax)/ MAX(1, MAX(fabs(attr1.Fmax), fabs(attr2.Fmax))));
		float fisim = 1 - fabs(attr1.Fint-attr2.Fint) /MAX(1, MAX(fabs(attr1.Fint), fabs(attr2.Fint)));
		//float lsim = 1 - fabs(attr1.Dlen-attr2.Dlen) / MAX(0.1, MAX(fabs(attr1.Dlen), fabs(attr2.Dlen)));
		float lsim = MIN(fabs(attr1.Dlen), fabs(attr2.Dlen)) / MAX(0.1f, MAX(fabs(attr1.Dlen), fabs(attr2.Dlen)));
		float asim = 1 - fabs(fabs(attr1.Dang) - fabs(attr2.Dang));
		float Q1minusQ2 = MAX(fabs(attr1.Q1-attr2.Q1), fabs(attr1.Q2-attr2.Q2));
		Q1minusQ2 = MAX(Q1minusQ2, fabs(attr1.Q3-attr2.Q3));
		float Q1norm = fabs(attr1.Q1 + attr1.Q2 + attr1.Q3) / 3;
		float Q2norm = fabs(attr2.Q1 + attr2.Q2 + attr2.Q3) / 3;
		//float qsim = MAX(-1, 1 - Q1minusQ2 / MAX(Q1norm, Q2norm));
		float qsim = 1- fabs(Q1norm - Q2norm) / MAX(1, MAX(Q1norm, Q2norm));
		qsim = (qsim < -1)? -1 : qsim;
		qsim = MIN(1, qsim);
		//attr1.print();
		//attr2.print();
		//printf("vsim = %f, tsim = %f, fisim = %f, isim = %f, lsim = %f, qsim = %f\n", vsim, tsim, fisim,
		//		isim, lsim, qsim);
		//return 0.25*vsim + 0.05*tsim + 0.49*fisim + 0.1*isim + 0.1*lsim + 0.01*qsim;
		return sim_w[0]*vsim + sim_w[1]*tsim + sim_w[2]*fisim + sim_w[3]*isim + sim_w[4]*lsim + sim_w[5]*qsim;
	} 	
	
	float simpleScore(DualNode* dualnode) {
		float vsim, fisim;
		if(MAX(volume, dualnode->volume) == 0) {
			vsim = 0;
		} else {
			vsim = MIN(volume, dualnode->volume) / MAX(volume, dualnode->volume);
		}
		float f1 = (volume == 0)? 0:(fint/volume);
		float f2 = (dualnode->volume == 0)? 0:(dualnode->fint/dualnode->volume);
		float len = MAX(1, MAX(fabs(f1), fabs(f2)));
		fisim = 1 - fabs(f2-f1)/len;
		return 0.5f*vsim + 0.5f*fisim;
	}
	

	/**
	 * stream output operator for DualNode.
	 */ 
	friend ostream& operator << (ostream& out, const DualNode& node) {
		int k;
		out << node.id << " " << node.level << endl;
		out << node.fmin << " " << node.fmax << " " << node.volume << " " << node.ncrit << " " << node.fint << endl;
		out << node.bl << node.bu << endl;
		out << node.moments << endl;
		out << node.lowerDegree() << " " << node.upperDegree() << " ";
		for(k = 0; k < node.lowerDegree(); k++) {
			out << (node.lowerNeighbor(k))->id << " ";
		}
		for(k = 0; k < node.upperDegree(); k++) {
			out << node.upperNeighbor(k)->id << " ";
		}
		out << endl;
		return out;
	}
	
	friend istream& operator >> (istream& in, DualNode& node) {
		in >> node.id >> node.level;
		in >> node.fmin >> node.fmax >> node.volume >> node.ncrit >> node.fint;
		in >> node.bl >> node.bu;
		in >> node.moments;
		return in;
	}
	/**
	 * Print the info about DualNode
	 */
	void print() {
		int k;
		printf("\t%d  %f, %f, [%f, %f], %d, (%d %d %d), (%d %d %d): ", id, volume, fint,
				fmin, fmax, ncrit, bl.b0, bl.b1, bl.b2, bu.b0, bu.b1, bu.b2);
		for (k = 0; k < lowerDegree(); k++) {
			printf(" %d ", lowerNeighbor(k)->id);
		}
		printf(" ---- ");
		for (k = 0; k < upperDegree(); k++) {
			printf(" %d ", upperNeighbor(k)->id);
		}
		printf("\n");
		//moments.printRaw();
#ifdef _MOMENTS 
		MomtAttrib att = moments.toAttributes();
		att.print();
#endif
	}
	
	void printPretty() {
		MomtAttrib att = moments.toAttributes();
		printf("%d\t %f\t %.2f, %.2f, %.2f, \t %.2f \t\t %.2f \t\t %.2f %.2f %.2f\n", id, volume, 
				att.I1, att.I2, att.I3, att.Fint, att.Dlen, att.Q1, att.Q2, att.Q3);
	}
	/*
	 *	Set the global min and max
	 */
	static void setGlobalMinMax(float & min, float & max) {
		all_min = min;
		all_max = max;
	}

	int id;
	int level;
	int ncrit;							// number of critical points in the volume
	float fmax, fmin;					// normalized [0, 1] max and min func values of this node
	float volume;						// normalized volume corresponding to this node;
	float fint;							// Integral of a second function 
	BettiNumber bl, bu;					// lower and upper boundary Betti numbers
	bool pruned;
#ifdef _MOMENTS	
	VolMoments moments;					// Volumetric moments
#endif

	std::vector<DualNode *> *p_below;
	std::vector<DualNode *> *p_above;
	DualNode* parent;					// points to the dual node in coarser level
	DualNode* match;					// A dual node in another graph that matches this node.		
private:
	static float all_min, all_max;		// global min and max
};

/**
 *	The class of Dual Contour Tree
 */
class DualGraph {
public:
	/**
	 * Construct a new dual contour tree
	 * @param nr: The number of functional ranges
	 * @param min: The minimum value of functional range
	 * @param max: Tha maximum value of functional range
	 */ 
	DualGraph(float min = 0, float max = 0, int nr = 0, float t_vol = 1.0f); 	
	
	//
	~DualGraph();

	/**
	 *	Build the dual graph from contour tree.
	 */
	void build(VolumeReg3Critical* p_data, AugmentedContourTree *actree);
	
	/**
	 * Load dual graph from a file
	 */
	bool loadDualFile(const char* fname);
	
	/**
	 *	Get the number of ranges
	 */
	int nRanges() const {
		return nrang;
	}

	/**
	 *	Construct a new dual graph by merging two adjacent subranges
	 */
	DualGraph* mergeRanges();

	/**
	 *	Print out info about the dual graph
	 */
	void print();

	/**
	 *	get the nth node in the level lvl
	 */
	DualNode* getNode(int lvl, int n) {
		return (nodes[lvl])[n];
	}

	/**
	 *	Get the total number of all nodes in the dual tree
	 */
	int totalNodes() const {
		return nnodes;
	}

	/**
	 * Get the number of nodes in the range i
	 */
	int nNodes(int i) const {
		return nodes[i].size();
	}

	/**
	 *	Save the dual graph into a file
	 */
	void dump(const char* fname);

	/**
	 *  Add a node to level i.
	 */
	void addDNote(int level, DualNode* p_dnode) {
#ifdef _DEBUG
		assert(level < nrang);
#endif
		nodes[level].insert(p_dnode);
	}

	/**
	 * Accessor to node j and level i.
	 */
	DualNode*& node(int level, int nid) {
#ifdef _DEBUG
		assert(level < nrang);
#endif		
		return nodes[level][nid];
	}

	/**
	 * Set the range of dual contour tree	
	 */
	void setRange(float min, float max) {
		m_min = min;
		m_max = max;
	}

	/**
	 * Remove nodes with volume < threshold
	 */
	void prune(float threshold);
	
	/**
	 * Insert a node to the level l of the tree
	 */
	void insert(int l, DualNode* pnode) {
		nodes[l].insert(pnode);
		nnodes ++;
	}

	float getTotalVolume() const {
		if(data == NULL) {
			fprintf(stderr, "\t\t getTotalVolume() data is NULL\n");
			return total_vol;
		}
		else return data->getVolume();
	}
	
	float getVolume() const { 
		return total_vol;
	}
	
	void setVolume(float t_vol) {
		total_vol = t_vol;
	}
	
protected:
	float	total_vol;
	float 	m_min, m_max;
	int		nrang;			// number of functional ranges
	int		nnodes;			// number of dual nodes
	float	*cut_vals;			// range cutting values
	dynamic_array<DualNode*> *nodes;	// Arrays of dual tree nodes
	dynamic_array<DualNode*> pruned;

	VolumeReg3Critical* data;

private:
	/**
	 * Release all memory
	 */
	void clean();
	
	/**
	 * Build a range of dual nodes
	 * @param ir: The range index
	 */
	void buildRange(int ir, dynamic_array<int>& crit_pnts, CutVertex* low_verts, int nl,
					CutVertex* up_verts, int nu, map<int, CutVertex*> *cut_maps, int levels[],
					VolumeReg3Critical* p_data, AugmentedContourTree* actree);

	/**
	 * get the cut level for a given value
	 */
	int getCutLevel(float val) {
		if(val < m_min) return -1;
		if(val > m_max) return nrang;
		int l = (int)MIN(nrang-1, (val-m_min)*nrang/(m_max-m_min));
		if(val <= cut_vals[l]) l --;
		return l;
	}
	 
};

#endif

