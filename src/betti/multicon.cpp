#include "multicon.h"
#include "pqueue.h"

#include <list>
#include <algorithm>

using namespace std;


//const static float VOL_THRESHOLD = 0.0001;
const static float VOL_THRESHOLD = 0.01f;


//#define _DEBUG 1
MultiConTree::MultiConTree(DualGraph *dtree, float total_vol)
{
	// we assume max_range is power of 2
	max_range = dtree->nRanges();
	max_level = log2(max_range);
	printf("range # = %d, level # = %d\n", max_range, max_level);
	dual_trees = (DualGraph **)malloc(sizeof(DualGraph *)*(max_level+1));
	dual_trees[max_level] = dtree;
	total_volume = total_vol;

	// create higher level dual graphs by merging ranges
	int i = max_level, j, k;
#ifdef _DEBUG
	dual_trees[i]->print();
#endif
	for (i = max_level; i > 0; i--) {
		dual_trees[i-1] = dual_trees[i]->mergeRanges();
#ifdef _DEBUG
		printf("dual tree at level %d\n", i-1);
		dual_trees[i-1]->print();
#endif
	}

	double vol = 0;
	for (i = 0; i < dual_trees[0]->nNodes(0); i++) {
		vol += dual_trees[0]->node(0, i)->volume;
	}
	fprintf(stdout, "Total volume: %f = (%f * %f)\n", vol*total_volume, vol, total_volume);
	for (i = 0; i <= max_level; i++) {
		for (j = 0; j < dual_trees[i]->nRanges(); j++) {
			for (k = 0; k < dual_trees[i]->nNodes(j); k++) {
				dual_trees[i]->node(j, k)->volume /= vol;
				dual_trees[i]->node(j, k)->fint /= vol;
			}
		}
	}

#ifdef _MOMENTS
/* 	float *eigval = fvector(1, 3);
	float **axes = matrix(1, 3, 1, 3);
	for (i = max_level; i >=0; i--) {
		for (j = 0; j < dual_trees[i]->nRanges(); j++) {
			printf("****************level %d****range %d************************\n", i, j);
			VolMoments mom;
			for (k = 0; k < dual_trees[i]->nNodes(j); k++) {
				mom = mom + dual_trees[i]->node(j, k)->moments;
			}
			mom.print();
			mom.principalAxes(eigval, axes);
			mom.toAttributes().print();
		}
	}
	free_vector(eigval, 1, 3);
	free_matrix(axes, 1, 3, 1, 3); */
#endif
	prune(VOL_THRESHOLD);

	for(int lev = max_level; lev >= 0; lev --) {
		char treeFile[64];
		sprintf(treeFile, "pruned_lev_%d.dct\0", lev);
		dual_trees[lev]->dump(treeFile);
	}
	//printf("level 0 dual tree after pruning\n");
	
#ifdef _COMPARE
	printf("ID \t Volume\t \t Moment of Inertia \t\t Pot Integral \t Dipole \t Moment of Quadrupole\n");
	for(i = 0; i <= 1; i++) {
		printf("----------------level %d------------------------\n", i);
		for(j = 0; j < dual_trees[i]->nRanges(); j++) {
			for(k = 0; k < dual_trees[i]->nNodes(j); k++) {
				if(dual_trees[i]->node(j, k)->volume >= 0.1*dual_trees[i]->nRanges()) {
					dual_trees[i]->node(j, k)->printPretty();
				}
			}
		}
	}
#else
	dual_trees[0]->print(); 
#endif

	FILE* log = fopen("DCT.log", "w");
	vector<float> volsize;
	for(i = 0; i < dual_trees[0]->nNodes(0) ; i++) {
		if(dual_trees[0]->node(0, i)->volume > 0.05f) 
			volsize.push_back(dual_trees[0]->node(0, i)->volume*vol*total_volume);
	}
	sort(volsize.begin(), volsize.end());
	for(i = 1, j = volsize.size(); j > 0 ; j--, i++) {
		fprintf(log, "%d %f\n", i, volsize[j-1]);
	}	
	fclose(log); 
	
}
void MultiConTree::prune(float threshold)
{
	dual_trees[0]->prune(threshold);

	// prune DCT at finer resoultions
	int level;
	for(level = 1; level <= max_level; level++) {
		dual_trees[level]->prune(threshold);
	}
}

MultiConTree::~MultiConTree()
{
	for (int i = 0; i <= max_level; i++) {
		delete dual_trees[i];
	}

	free(dual_trees); 
}

void MultiConTree::getOrientation(int l, int r, float *ctr, float **axes)
{
#ifdef _MOMENTS
	assert(l <= max_level && r < dual_trees[l]->nRanges());

	float *eigval = fvector(1, 3);
	float **eigvec = matrix(1, 3, 1, 3);
	VolMoments mom;
	for (int k = 0; k < dual_trees[l]->nNodes(r); k++) {
		mom = mom + dual_trees[l]->node(r, k)->moments;
	}	
	mom.print();
	mom.getCOM(ctr);
	mom.principalAxes(eigval, eigvec);
	for(int i = 0; i < 3; i++) {
		for(int j = 0; j < 3; j++) {
			axes[i][j] = eigvec[i+1][j+1];
		}
	}
	free_vector(eigval, 1, 3);
	free_matrix(eigvec, 1, 3, 1, 3);
#endif
}

float MultiConTree::matching(MultiConTree* mulcon)
{
	int level = MIN(maxLevel(), mulcon->maxLevel());
	list<DualNode *>* matched_lists = new list<DualNode *>[level+1];


	int i, j, k;
	PriorityQueue<DualNode *> *que;
	DualNode *dualnode;
	float tvol[2];

	tvol[0] = totalVolume();
	tvol[1] = mulcon->totalVolume();
#ifdef _DEBUG
	printf("total volume of tree 1 = %f, of tree 2 = %f\n", tvol[0], tvol[1]);
#endif

	que = new PriorityQueue<DualNode *>;
	// First try to match the level 0
	for (i = 0; i < dual_trees[0]->nNodes(0); i++) {
		dualnode = dual_trees[0]->node(0, i);
		que->insert(dualnode, dualnode->weight());
	}
	while (!que->isEmpty()) {
		que->extract(dualnode);
		DualNode *mtch_node = matchSearch(dualnode, mulcon->dual_trees[0]);
		if (mtch_node != NULL) {
			dualnode->setMatchNode(mtch_node);
			matched_lists[0].push_back(dualnode);
		}
	}
	//(dual_trees[0]->getNode(0, 0))->setMatchNode(mulcon->dual_trees[0]->getNode(0, 0));
	//matched_lists[0].push_back(dual_trees[0]->getNode(0, 0));


	for (i = 1; i <= level; i++) {
		for (j = 0; j < dual_trees[i]->nRanges(); j++) {
			for (k = 0; k < dual_trees[i]->nNodes(j); k++) {
				dualnode = dual_trees[i]->node(j, k);
				if (dualnode->parent->isMatched())
					que->insert(dualnode, dualnode->weight());
			}
		}
		while (!que->isEmpty()) {
			que->extract(dualnode); 
			DualNode *mtch_node = matchSearch(dualnode, mulcon->dual_trees[i]);
			//assert(mtch_node != NULL);
			if (mtch_node != NULL) {
				dualnode->setMatchNode(mtch_node);
				matched_lists[i].push_back(dualnode);
			}
		}
	}
	delete que;

	float score = 0;
	// compute the similarity score between two multi-res con trees
	for (i = 0; i <= level; i++) {
		//printf("---------------------------------level %d scores---------------------------\n", i);
		while (!matched_lists[i].empty()) {
			dualnode = matched_lists[i].front();
			matched_lists[i].pop_front();
			float s = dualnode->simScore2(dualnode->match);
			//printf("match volume: %f -> %f\n", dualnode->volume, dualnode->match->volume);
			//dualnode->print();
			//dualnode->match->print();
			//printf("s = %f\n", s);
			score += s * (dualnode->volume + dualnode->match->volume) / 2;
		}
	}
	score /= (level+1)*MAX(tvol[0], tvol[1]);
	//score /= level*MAX(tvol[0], tvol[1]);
	delete[] matched_lists;

	return score;
}

/************************************************************************/
/* Private functions                                                    */
/************************************************************************/
int MultiConTree::log2(int n) {
	assert(n > 0);
	int k = 0;
	while (n > 1) {
		n >>= 1 ;
		k++;
	}
	return k;
}

DualNode* MultiConTree::matchSearch(DualNode* dualnode, DualGraph* graph)
{
	int i, lvl = dualnode->level;
	PriorityQueue<DualNode *> candidates;
	DualNode* mtch_node = NULL;
	float grade;
	for (i = 0; i < graph->nNodes(lvl); i++) {
		DualNode *p_node = graph->node(lvl, i);
		if (p_node->isParentMatched(dualnode) && !(p_node->isMatched())) {
			candidates.insert(graph->node(lvl, i), graph->node(lvl, i)->weight());
			//printf("*************candidtate: %d %f\n", (graph->nodes[lvl])[i]->id, (graph->nodes[lvl])[i]->weight());
		}
	}
	while (!candidates.isEmpty()) {
		DualNode *cand;
		candidates.extract(cand);
		float scor = matchScore(dualnode, cand);
		if (mtch_node == NULL) {
			grade = scor;
			mtch_node = cand;
		} else if (scor > grade) {
			grade = scor;
			mtch_node = cand;
		}
	}
	return mtch_node;
}

float MultiConTree::matchScore(DualNode *node1, DualNode *node2)
{
	float s;
	//s = node1->simScore(node2) - (node1->simScore(node1) + node2->simScore(node2)) / 2;
	//s = node1->simScore2(node2) - (node1->simScore2(node1) + node2->simScore2(node2)) / 2;
	s = node1->simpleScore(node2);
	return s;
}
