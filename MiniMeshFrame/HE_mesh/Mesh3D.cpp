#include "Mesh3D.h"

#include <fstream>
#include <iostream>
#include <xutility>

#define SWAP(a,b,T) {T tmp=(a); (a)=(b); (b)=tmp;}
#define min(a,b) a<b?a:b
#define max(a,b) a>b?a:b


Mesh3D::Mesh3D(void)
{
	// intialization
	pvertices_list_ = NULL;
	pfaces_list_ = NULL;
	pedges_list_ = NULL;

	xmax_ = ymax_ = zmax_ = 1.f;
	xmin_ = ymin_ = zmin_ = -1.f;

	num_components_ = 0;
	average_edge_length_ = 1.f;
}

void Mesh3D::ClearData(void)
{
	ClearVertex();
	ClearEdges();
	ClearFaces();
	edgemap_.clear();

	xmax_ = ymax_ = zmax_ = 1.f;
	xmin_ = ymin_ = zmin_ = -1.f;
}

void Mesh3D::ClearVertex(void)
{

	if (pvertices_list_==NULL)
	{
		return;
	}
	else
	{
		for (VERTEX_ITER viter = pvertices_list_->begin(); viter != pvertices_list_->end(); viter++)
		{
			if (*viter != NULL)
			{
				delete *viter;
				*viter = NULL;
			}
			else
			{
				// ERROR
			}
		}
		delete pvertices_list_;
		pvertices_list_ = NULL;
	}
}

void Mesh3D::ClearEdges(void)
{
	if (pedges_list_ == NULL)
	{
		return;
	}
	else
	{
		for (EDGE_ITER eiter = pedges_list_->begin(); eiter!=pedges_list_->end(); eiter++)
		{
			if (*eiter != NULL)
			{
				delete *eiter;
				*eiter = NULL;
			}
			else
			{
				// ERROR
			}
		}
		delete pedges_list_;
		pedges_list_ = NULL;
	}
}

void Mesh3D::ClearFaces(void)
{
	if (pfaces_list_==NULL)
	{
		return;
	}
	else
	{
		for (FACE_ITER fiter = pfaces_list_->begin(); fiter!=pfaces_list_->end(); fiter++)
		{
			if (*fiter != NULL)
			{
				delete *fiter;
				*fiter = NULL;
			}
			else
			{
				// ERROR
			}
		}
		delete pfaces_list_;
		pfaces_list_ = NULL;
	}
}

HE_vert* Mesh3D::InsertVertex(const Vec3f& v)
{
	HE_vert* pvert = new HE_vert(v);
	if (pvertices_list_ == NULL)
	{
		pvertices_list_ = new std::vector<HE_vert*>;
	}
	pvert->id_ = static_cast<int>(pvertices_list_->size());
	pvertices_list_->push_back(pvert);
	return pvert;
}

HE_edge* Mesh3D::InsertEdge(HE_vert* vstart, HE_vert* vend)
{
	if (vstart==NULL || vend==NULL)
	{
		return NULL;
	}

	if (pedges_list_==NULL)
	{
		pedges_list_ = new std::vector<HE_edge*>;
	}

	if (edgemap_[PAIR_VERTEX(vstart, vend)] != NULL)
	{
		return edgemap_[PAIR_VERTEX(vstart, vend)];
	}

	HE_edge* pedge = new HE_edge;
	pedge->pvert_ = vend;
	pedge->start_ = vstart;
	//pedge->pvert_->degree_ ++;
	vstart->pedge_ = pedge;
	//edgemap_[PAIR_VERTEX(vstart, vend)] = pedge;

	pedge->id_ = static_cast<int>(pedges_list_->size());
	pedges_list_->push_back(pedge);

	return pedge;
}

HE_face* Mesh3D::InsertFace(std::vector<HE_vert* >& vec_hv)
{
	int vsize = static_cast<int>(vec_hv.size());
	//if (vsize != 3)
	//{
	//	return NULL;
	//}

	if (pfaces_list_ == NULL)
	{
		pfaces_list_ = new std::vector<HE_face*>;
	}

	HE_face *pface = new HE_face;
	pface->valence_ = vsize;
	VERTEX_ITER viter = vec_hv.begin();
	VERTEX_ITER nviter = vec_hv.begin();
	nviter ++;

	HE_edge *he1=NULL, *he2=NULL;
	std::vector<HE_edge*> vec_edges;
	int i=0;
	for (i=0; i<vsize-1; i++)
	{
		he1 = InsertEdge( *viter, *nviter);
		he2 = InsertEdge( *nviter, *viter);

		if (pface->pedge_==NULL) 
			pface->pedge_ = he1;

		he1->pface_ = pface;
		he1->ppair_ = he2;
		he2->ppair_ = he1;
		vec_edges.push_back(he1);
		viter++, nviter++;
	}

	nviter = vec_hv.begin();

	he1 = InsertEdge(*viter, *nviter);
	he2 = InsertEdge(*nviter , *viter);
	he1->pface_ = pface;
	if (pface->pedge_==NULL) 
		pface->pedge_ = he1;

	he1->ppair_ = he2;
	he2->ppair_ = he1;
	vec_edges.push_back(he1);

	for (i=0; i<vsize-1; i++) 
	{
		vec_edges[i]->pnext_ = vec_edges[i+1];
		vec_edges[i+1]->pprev_ = vec_edges[i];
	}
	vec_edges[i]->pnext_ = vec_edges[0];
	vec_edges[0]->pprev_ = vec_edges[i];

	pface->id_ = static_cast<int>(pfaces_list_->size());
	pfaces_list_->push_back(pface);

	return pface;
}

bool Mesh3D::LoadFromOBJFile(const char* fins)
{

	if (fins == NULL)
	{
		return false;
	}
	FILE *pfile = fopen(fins, "r");
	fseek(pfile, 0, SEEK_SET);
	char pLine[512];
// 	//if (pvertices_list_ == NULL)
// 	{
// 		//pvertices_list_ = new std::vector<HE_vert>;
// 		//printf("num of thread,%d,list address,%d,\r\n",omp_get_thread_num(),pvertices_list_);
// 	}
	while (fgets(pLine, 512, pfile))
	{
		char *p[10] = { NULL };
		char *saveptr = NULL;
		p[0] = strtok_s(pLine, " ", &saveptr);
		//HE_vert* cur_ = NULL;
		Vec3f* nvv = new Vec3f;
		(*nvv)[0] = (double)atof(p[0]);
		for (int i = 1; p[i] = strtok_s(NULL, " ", &saveptr); i++)
		{

			if (i < 2)
			{
				(*nvv)[i] = (double)atof(p[i]);
			}
			else if (i == 2)
			{
				(*nvv)[i] = (double)atof(p[i]);
				InsertVertex((*nvv));
 				(*(pvertices_list_->rbegin()))->set_seleted(UNSELECTED);// all point init with 0;
				(*(pvertices_list_->rbegin()))->sumVector = (0.0, 0.0, 0.0);
				(*(pvertices_list_->rbegin()))->truth = true;
				(*(pvertices_list_->rbegin()))->subgraphflag=-1;
			}
			else if (*p[i] != '\n')
			{
				(*(pvertices_list_->rbegin()))->neighborIdx.push_back((int)atof(p[i]) - 1);
				(*(pvertices_list_->rbegin()))->neighbor_search_.push_back((int)atof(p[i]) - 1);
				(*(pvertices_list_->rbegin()))->degree_++;
			}
		}
		delete nvv;
	}
	fclose(pfile);
	//return true;
	//num_of_vertex_list() = pvertices_list_->size();
	for (auto iter = pvertices_list_->begin(); iter != pvertices_list_->end(); iter++)
	{
		if ((*iter)->degree() == 1)
		{
			(*iter)->set_boundary_flag(BOUNDARY);
			findLooppoint((*iter)->id());
		}
		else
		{
			(*iter)->set_boundary_flag(INNER);
		}
		

	}
	for (auto iter = pvertices_list_->begin(); iter != pvertices_list_->end(); iter++)
	{
		if ((*iter)->selected() == UNSELECTED)
		{
			(*iter)->set_boundary_flag(IN_LOOP);
		}
		else
			(*iter)->set_seleted(UNSELECTED);
	}
	num_vertex_ = pvertices_list_->size();
	return true;
}

void Mesh3D::WriteToOBJFile(const char* fouts)
{
	std::ofstream fout(fouts);

	fout<<"g object\n";
	fout.precision(16);
	//output coordinates of each vertex
	VERTEX_ITER viter = pvertices_list_->begin();
	for (;viter!=pvertices_list_->end(); viter++) 
	{
		fout<<"v "<< std::scientific <<(*viter)->position_.x() 
			<<" "<<(*viter)->position_.y() <<" "<< (*viter)->position_.z() <<"\n";
	}

	// 		for (viter = pvertices_list_->begin();viter!=pvertices_list_->end(); viter++) 
	// 		{
	// 			fout<<"vn "<< std::scientific <<(*viter)->normal_.x() 
	// 				<<" "<<(*viter)->normal_.y() <<" "<<(*viter)->normal_.z() <<"\n";
	// 		}
	//output the valence of each face and its vertices_list' id

	FACE_ITER fiter = pfaces_list_->begin();

	for (;fiter!=pfaces_list_->end(); fiter++) 
	{
		fout<<"f";

		HE_edge* edge = (*fiter)->pedge_; 

		do {
			fout<<" "<<edge->ppair_->pvert_->id_+1;
			edge = edge->pnext_;

		} while (edge != (*fiter)->pedge_);
		fout<<"\n";
	}

	fout.close();
}

void Mesh3D::UpdateMesh(void)
{
	if (!isValid())
	{
		std::cout << "Invalid" << "\n";
		return;
	}
	SetBoundaryFlag();
	BoundaryCheck();
	UpdateNormal();
	ComputeBoundingBox();
	ComputeAvarageEdgeLength();
	SetNeighbors();
}

void Mesh3D::SetBoundaryFlag(void)
{
	for (EDGE_ITER eiter = pedges_list_->begin(); eiter!=pedges_list_->end(); eiter++)
	{
		if ((*eiter)->pface_ == NULL)
		{
			(*eiter)->set_boundary_flag(BOUNDARY);
			(*eiter)->ppair_->set_boundary_flag(BOUNDARY);
			(*eiter)->pvert_->set_boundary_flag(BOUNDARY);
			(*eiter)->ppair_->pvert_->set_boundary_flag(BOUNDARY);
			(*eiter)->ppair_->pface_->set_boundary_flag(BOUNDARY);
		}
	}
}

void Mesh3D::BoundaryCheck()
{
	for (VERTEX_ITER viter=pvertices_list_->begin(); viter!=pvertices_list_->end(); viter++)
	{
		if ((*viter)->isOnBoundary())
		{
			HE_edge* edge = (*viter)->pedge_;
			int deg = 0;
			while (edge->pface_!=NULL && deg<(*viter)->degree())
			{
				edge = edge->pprev_->ppair_;
				deg ++;
			}
			(*viter)->pedge_ = edge;
		}
	}
}

void Mesh3D::UpdateNormal(void)
{
	ComputeFaceslistNormal();
	ComputeVertexlistNormal();
}

void Mesh3D::ComputeFaceslistNormal(void)
{
	for (FACE_ITER fiter = pfaces_list_->begin(); fiter!=pfaces_list_->end(); fiter++)
	{
		ComputePerFaceNormal(*fiter);
	}
}

void Mesh3D::ComputePerFaceNormal(HE_face* hf)
{
	HE_edge *pedge = hf->pedge_;
	HE_edge *nedge = hf->pedge_->pnext_;

	HE_vert *p = pedge->pvert_;
	HE_vert *c = pedge->pnext_->pvert_;
	HE_vert *n = nedge->pnext_->pvert_;

	Vec3f pc, nc;
	pc = p->position_ - c->position_;
	nc = n->position_ - c->position_;

	hf->normal_ = nc ^ pc;	// cross prodoct
	hf->normal_.normalize();
}

void Mesh3D::ComputeVertexlistNormal(void)
{
	for (VERTEX_ITER viter = pvertices_list_->begin(); viter!=pvertices_list_->end(); viter++) 
	{
		ComputePerVertexNormal(*viter);
	}
}

void Mesh3D::ComputePerVertexNormal(HE_vert* hv)
{
	if (hv->degree_ < 2)
	{
		// ERROR: the degree of the vertex is less than 2
		hv->normal_ = Vec3f(1.f,0.f,0.f);
		return;
	}

	HE_edge *edge = hv->pedge_;
	if (edge == NULL)
	{
		// ERROR: the edge attached to the vertex is NULL
		hv->normal_ = Vec3f(1.f,0.f,0.f);
		return;
	}

	hv->normal_ = Vec3f(0.f,0.f,0.f);
	if (hv->boundary_flag_ == INNER)
	{
		int iterNum = 0;
		do 
		{
			iterNum++;
			if (iterNum > hv->degree())
			{
				/*hv->set_position(hv->position() * 1.1f);*/
				std::cout << "    iterNum > hv->degree : " << hv->id() << "\n";
				break;
			}
			//hv->normal_ = hv->normal_ + edge->pface_->normal_;
			Vec3f  p = edge->pvert_->position(),
				q = edge->pnext_->pvert_->position(),
				r = edge->pprev_->pvert_->position();
			Vec3f  n = (q-p) ^ (r-p);
			hv->normal_ = hv->normal_ + n;
			edge = edge->ppair_->pnext_;
		} while (edge != hv->pedge_ && edge != NULL);
	}
	else
	{
		// NOTE: for the boundary vertices, this part may be something wrong
		//	     Up to now, define the normals all as unity
		hv->normal_ = Vec3f(1.f, 0.f, 0.f);

		//int degree_flag = 0;
		//for (int i=0; i<hv->degree_-1; i++)
		//{
		//	edge = edge->ppair_->pnext_;
		//	if (edge == NULL)
		//	{
		//		// ERROR: the algorithm of computing boundary vertices has errors!
		//		break;
		//	}
		//	if (edge->pface_ != NULL)
		//	{
		//		hv->normal_ = hv->normal_ + edge->pface_->normal_;
		//	}
		//}
	}
	hv->normal_.normalize();
}

void Mesh3D::ComputeBoundingBox(void)
{
	if (pvertices_list_->size() < 3)
	{
		return;
	}

#define MAX_FLOAT_VALUE (static_cast<float>(10e10))
#define MIN_FLOAT_VALUE	(static_cast<float>(-10e10))
	
	xmax_ = ymax_ = zmax_ = MIN_FLOAT_VALUE;
	xmin_ = ymin_ = zmin_ = MAX_FLOAT_VALUE;

	VERTEX_ITER viter = pvertices_list_->begin();
	for (; viter!=pvertices_list_->end(); viter++)
	{
		xmin_ = min(xmin_, (*viter)->position_.x());
		ymin_ = min(ymin_, (*viter)->position_.y());
		zmin_ = min(zmin_, (*viter)->position_.z());
		xmax_ = max(xmax_, (*viter)->position_.x());
		ymax_ = max(ymax_, (*viter)->position_.y());
		zmax_ = max(zmax_, (*viter)->position_.z());
	}
}

void Mesh3D::Unify(float size)
{
	float scaleX = xmax_ - xmin_;
	float scaleY = ymax_ - ymin_;
	float scaleZ = zmax_ - zmin_;
	float scaleMax;

	if (scaleX < scaleY)
	{
		scaleMax = scaleY;
	}
	else
	{
		scaleMax = scaleX;
	}
	if (scaleMax < scaleZ)
	{
		scaleMax = scaleZ;
	}
	float scaleV = size / scaleMax;
	Vec3f centerPos((xmin_ + xmax_) / 2.f, (ymin_ + ymax_) / 2.f, (zmin_ + zmax_) / 2.f);
	for (size_t i = 0; i != pvertices_list_->size(); i++)
	{
		pvertices_list_->at(i)->position_ = (pvertices_list_->at(i)->position_ - centerPos) * scaleV;
	}
}

void Mesh3D::ComputeAvarageEdgeLength(void)
{
	if(!isValid())
	{
		average_edge_length_ = 0.f;
		return;
	}
	float aveEdgeLength = 0.f;
	for (int i=0; i<num_of_half_edges_list(); i++)
	{
		HE_edge* edge = get_edges_list()->at(i);
		HE_vert* v0 = edge->pvert_;
		HE_vert* v1 = edge->ppair_->pvert_;
		aveEdgeLength += (v0->position() - v1->position()).length();
	}
	average_edge_length_ = aveEdgeLength/num_of_half_edges_list();
	//std::cout << "Average_edge_length = " << average_edge_length_ << "\n";
}

HE_face* Mesh3D::get_face(int vId0, int vId1, int vId2)
{
	HE_vert *v0 = get_vertex(vId0);
	HE_vert *v1 = get_vertex(vId1);
	HE_vert *v2 = get_vertex(vId2);
	if (!v0 || !v1 || !v2)
	{
		return NULL;
	}

	HE_face* face=NULL;

	// 由于对边界点的邻域遍历有bug，所以找到非边界点进行邻域遍历
	if (v0->isOnBoundary())
	{
		if (!v1->isOnBoundary())
		{
			SWAP(v0, v1, HE_vert*);
		}
		else if (!v2->isOnBoundary())
		{
			SWAP(v0, v2, HE_vert*);
		}
		else
		{
			// v0, v1, v2 都是边界点
			// 暂时先不处理
			return NULL;
		}
	}

	if (!v0->isOnBoundary())	// 对边界点的遍历有bug
	{
		HE_edge* edge=v0->pedge_;
		bool inFace = true;
		do 
		{
			bool b1 = isFaceContainVertex(edge->pface_, v1);
			bool b2 = isFaceContainVertex(edge->pface_, v2);
			if (!b1 && !b1)
			{
				edge = edge->ppair_->pnext_;
			}
			else if(b1 && b2)
			{
				face = edge->pface_;
				break;
			}
			else
			{
				inFace = false;
				break;
			}
		} while (edge!=v0->pedge_ && edge!=NULL);
	}

	return face;
}

HE_face* Mesh3D::get_face(const std::vector<unsigned int>& ids)
{
	if (ids.size()<3)
	{
		std::cout << "查询点数过少，无法返回面\n";
		return NULL;
	}
	// 首先找到一个非边界点
	HE_vert* v = NULL;
	for (unsigned int i=0; i<ids.size(); i++)
	{
		if (!get_vertex(ids[i])->isOnBoundary())
		{
			v = get_vertex(ids[i]);
			break;
		}
	}
	if (!v)
	{
		// 所有点都是边界点
		// 暂不处理
		return NULL;
	}

	HE_edge *edge = v->pedge_;
	HE_face *face = NULL;
	do 
	{
		face = edge->pface_;
		edge = edge->ppair_->pnext_;
		bool bInFace = isFaceContainVertex(face, get_vertex(ids[0]));
		if (!bInFace)
		{
			continue;
		}
		for (unsigned int i=1; i<ids.size(); i++)
		{
			bool b = isFaceContainVertex(face, get_vertex(ids[i]));
			if (b!=bInFace)
			{
				bInFace = false;
				break;
			}
		}
		if (bInFace)
		{
			return face;
		}
	} while (edge!=v->pedge_ && edge!=NULL);
	return NULL;
}

bool Mesh3D::isFaceContainVertex(HE_face* face, HE_vert* vert)
{
	HE_edge* edge = face->pedge_;
	do 
	{
		if (edge->pvert_==vert)
		{
			return true;
		}
		edge = edge->pnext_;
	} while (edge!=face->pedge_ && edge!=NULL);
	return false;
}

int Mesh3D::GetFaceId(HE_face* face)
{
	return !face ? -1 : face->id();
}

void Mesh3D::ResetFaceSelectedTags(int tag)
{
	for (int i=0; i<num_of_face_list(); i++)
	{
		get_face(i)->set_selected(tag);
	}
}

void Mesh3D::ResetVertexSelectedTags(int tag)
{
	for (int i=0; i<num_of_vertex_list(); i++)
	{
		get_vertex(i)->set_seleted(tag);
	}
}

bool Mesh3D::isNeighbors(HE_vert* v0, HE_vert* v1)
{
	if (!v0 || !v1)
	{
		return false;
	}

	HE_edge *edge = v0->pedge_;
	do 
	{
		if (edge->pvert_==v1)
		{
			return true;
		}
		edge = edge->ppair_->pnext_;
	} while (edge!=v0->pedge_ && edge);
	return false;
}

int Mesh3D::GetSelectedVrtId()
{
	if (!isValid())
	{
		return -1;
	}
	for (int i=0; i<num_of_vertex_list(); i++)
	{
		if (get_vertex(i)->selected()==SELECTED)
		{
			return i;
		}
	}
	return -1;
}

void Mesh3D::CreateMesh(const std::vector<Vec3f>& verts, const std::vector<int>& triIdx)
{
	ClearData();
	for (unsigned int i=0; i<verts.size(); i++)
	{
		InsertVertex(verts[i]);
	}
	for (unsigned int i=0; i<triIdx.size(); i=i+3)
	{
		std::vector<HE_vert*> tri;
		HE_vert *v0 = get_vertex(triIdx[i]);
		HE_vert *v1 = get_vertex(triIdx[i+1]);
		HE_vert *v2 = get_vertex(triIdx[i+2]);
		if (!v0 || !v1 || !v2) continue;
		tri.push_back(v0);
		tri.push_back(v1);
		tri.push_back(v2);
		InsertFace(tri);
	}
	UpdateMesh();
}

void Mesh3D::CreateMesh(const std::vector<double>& verts, const std::vector<unsigned>& triIdx)
{
	ClearData();
	for (unsigned int i=0; i<verts.size(); i=i+3)
	{
		InsertVertex(Vec3f(verts[i], verts[i+1], verts[i+2]));
	}
	for (unsigned int i=0; i<triIdx.size(); i=i+3)
	{
		std::vector<HE_vert*> tri;
		HE_vert *v0 = get_vertex(triIdx[i]);
		HE_vert *v1 = get_vertex(triIdx[i+1]);
		HE_vert *v2 = get_vertex(triIdx[i+2]);
		if (!v0 || !v1 || !v2) continue;
		tri.push_back(v0);
		tri.push_back(v1);
		tri.push_back(v2);
		InsertFace(tri);
	}
	UpdateMesh();
}

int Mesh3D::GetBoundaryVrtSize()
{
	int count = 0;
	for (int i=0; i<num_of_vertex_list(); i++)
	{
		if (get_vertex(i)->isOnBoundary())
		{
			count ++;
		}
	}
	return count;
}

Mesh3D::~Mesh3D(void)
{
	ClearData();
}
bool HE_vert::set_cutchoice(int cutchoiceID)
{

	for (auto iter = unavailiable.begin(); iter != unavailiable.end(); iter++)
	{
		if (cutchoiceID == *iter)
			return false;
	}
	cutchoice.push_back(cutchoiceID);
	last_unavailiable_size.push_back(unavailiable.size());
	for (int t1 = 0; t1 < cutedpoints.size(); t1++)//¶ÔÓÚÃ¿ÖÖÇÐžî·œÊœ
	{
		for (int t2 = 0; t2 < cutedpoints.at(t1).size(); t2++)//¶ÔÓÚÇÐžî·œÊœÖÐµÄÃ¿Ò»Ìõ±ß
		{
			for (int t3 = 0; t3 < cutedpoints.at(cutchoiceID).size(); t3++)//¶ÔÓÚÑ¡ÔñµÄÇÐžî·œÊœµÄÃ¿Ò»Ìõ±ß
			{
				if (cutedpoints.at(t1).at(t2) == cutedpoints.at(cutchoiceID).at(t3))//Èç¹ûÒÑŸ­ÓÐÏàÍ¬µÄ±ßÁË
				{
					//	printf("id %d,unavailable size %d,thread%d\n",this->id(), unavailiable.size(),omp_get_thread_num());
					unavailiable.push_back(t1);
					break;
				}
			}
		}
	}
	//printf("*********************************************\n");
	return true;
};
bool HE_vert::reset_cutchoice() {
	for (int w = unavailiable.size(); w > last_unavailiable_size.back(); w--)
	{
		unavailiable.pop_back();
	}
	last_unavailiable_size.pop_back();
	return true;
}

void Mesh3D::findLooppoint(int ID)
{
	int Degreee = pvertices_list_->at(ID)->degree();
	for (auto iter1 = pvertices_list_->at(ID)->neighborIdx.begin(); iter1 != pvertices_list_->at(ID)->neighborIdx.end(); iter1++)
	{
		if (pvertices_list_->at(*iter1)->selected() == SELECTED)
		{
			//cout << "the neighbor "<<*iter1<< " has been set SELECTED"<<endl;
			--Degreee;
			//cout << Degreee<< endl;
		}
	}
	if (Degreee <= 1)
	{
		pvertices_list_->at(ID)->set_seleted(SELECTED);
		//cout << ID << endl;
		for (auto iter2 = pvertices_list_->at(ID)->neighborIdx.begin(); iter2 != pvertices_list_->at(ID)->neighborIdx.end(); iter2++)
		{
			if (pvertices_list_->at(*iter2)->selected() == UNSELECTED)
			{
				findLooppoint(pvertices_list_->at(*iter2)->id());
			}
		}
	}
	else
	{
		return;
	}
}

bool Mesh3D::computeRoot()
{
	std::ofstream outfile("result.txt");
	if (!outfile) 
	{ 
		std::cout << "file open wrong\n"; return false; 
	}

	int truth = 0;
	//#pragma omp parallel for shared(truth) 
	for (int loop_0_ = 0; loop_0_ < pvertices_list_->size(); loop_0_++)
	{
		if (pvertices_list_->at(loop_0_)->boundary_flag() == BOUNDARY)continue;
		for (int loop_1_ = 0; loop_1_ < pvertices_list_->size(); loop_1_++)
		{
			if (pvertices_list_->at(loop_1_)->boundary_flag() == BOUNDARY)continue;
			for (int loop_2_ = 0; loop_2_ < pvertices_list_->size(); loop_2_++)
			{
				if (pvertices_list_->at(loop_2_)->boundary_flag() == BOUNDARY)continue;
				for (int loop_3_ = 0; loop_3_ < pvertices_list_->size(); loop_3_++)
				{
					if (pvertices_list_->at(loop_3_)->boundary_flag() == BOUNDARY)continue;

					Computecutchoice(loop_0_);
					Computecutchoice(loop_1_);
					Computecutchoice(loop_2_);
					Computecutchoice(loop_3_);
					int root[4] = { loop_0_,loop_1_,loop_2_,loop_3_ };
					//cout<<root[0]<<" "<<root[1]<<" "<<root[2]<<" "<<root[3]<<endl;
					//continue;
					//if (omp_get_thread_num() == 0)
						//outfile << root[0] << " " << root[1] << " " << root[2] << " " << root[3] << "\n";
					//printf("num of thread,%d,\r\n",omp_get_num_threads());
					//printf("num of thread,%d,%d,%d,%d,%d,\r\n", omp_get_thread_num(), root[0], root[1], root[2], root[3]);
					for (int i = 0; i < 4; i++)
					{
						pvertices_list_->at(root[i])->cutchoice.clear();
						pvertices_list_->at(root[i])->unavailiable.clear();
						pvertices_list_->at(root[i])->last_unavailiable_size.clear();
					}

					for (int num_cut_0_ = 0; num_cut_0_ < pvertices_list_->at(loop_0_)->cutedpoints.size(); num_cut_0_++)
					{
						if (pvertices_list_->at(root[0])->set_cutchoice(num_cut_0_) == false)continue;

						for (int num_cut_1_ = 0; num_cut_1_ < pvertices_list_->at(loop_1_)->cutedpoints.size(); num_cut_1_++)
						{
							if (pvertices_list_->at(root[1])->set_cutchoice(num_cut_1_) == false)continue;

							for (int num_cut_2_ = 0; num_cut_2_ < pvertices_list_->at(loop_2_)->cutedpoints.size(); num_cut_2_++)
							{
								if (pvertices_list_->at(root[2])->set_cutchoice(num_cut_2_) == false)continue;

								for (int num_cut_3_ = 0; num_cut_3_ < pvertices_list_->at(loop_3_)->cutedpoints.size(); num_cut_3_++)
								{
									if (pvertices_list_->at(root[3])->set_cutchoice(num_cut_3_) == false)continue;


									int numcuts[4] = { num_cut_0_,num_cut_1_,num_cut_2_,num_cut_3_ };
									if (!findsubgraph(root, 4, numcuts, 4))
									{
										Reset(root,4);//reset pvertices list
										pvertices_list_->at(root[3])->reset_cutchoice();
										continue;
									}
									
									//return true;
									if (computeSubgraph(0))
									{
										if (computeSubgraph(1))
										{
											if (computeSubgraph(2))
											{
												if (computeSubgraph(3))
												{
//#pragma omp critical(second)
													{
														//printf("root %d,root %d,root %d,root %d,thread %d", root[0], root[1], root[2], root[3], omp_get_thread_num());
														//outfile << omp_get_num_threads() << "\n";
														//outfile << omp_get_thread_num() << "\n";
														outfile << "split skeleton successfully!" << "\n";
														outfile << root[0] << " " << root[1] << " " << root[2] << " " << root[3] << "\n";
														outfile << numcuts[0] << " " << numcuts[1] << " " << numcuts[2] << " " << numcuts[3] << "\n";
														outfile.close();
//#pragma omp atomic
														truth++;
													}
												}
											}
										}
									}
									if (truth)break;
									Reset(root,4);
									pvertices_list_->at(root[3])->reset_cutchoice();
								}
								pvertices_list_->at(root[2])->reset_cutchoice();
								if (truth)break;
							}
							pvertices_list_->at(root[1])->reset_cutchoice();
							if (truth)break;
						}
						pvertices_list_->at(root[0])->reset_cutchoice();
						if (truth)break;
					}

					if (truth)break;
				}
				if (truth)	break;
			}
			if (truth)break;
		}
		if (truth)loop_0_ = pvertices_list_->size();
	}
	if (!truth)
	{
		outfile << "split failed, cannot find 4 root to fit the needs" << std::endl;
		outfile.close();
	}
	return false;
}

int Mesh3D::Computecutchoice(int i)
{

	HE_vert& root = *(pvertices_list_->at(i));
	//root.cutedpoints.clear();//for the usecompute
	if (root.cutedpoints.size() > 0 || root.neighborIdx.size() == 0)
	{
		return root.cutedpoints.size();
	}

	for (int k = 0; k < root.neighborIdx.size() - 1; k++)
	{
		std::vector<int> cutedpoint;
		if (k == 0)
		{
			for (unsigned int j = 0; j < root.neighborIdx.size(); j++)
			{
				cutedpoint.clear();
				//cutedpoint.push_back(root.neighbor_search_.at(j));
				cutedpoint.push_back(j);//store the index of neighbor array
				root.cutedpoints.push_back(cutedpoint);
			}
		}
		else
		{
			std::vector<std::vector<int>> p = root.cutedpoints;
			for (auto iter = p.begin(); iter != p.end(); iter++)//ѡ\D4\F1һ\D6\D6\D2Ѿ\AD\BC\C6\CB\E3\BAõ\C4cut·\BE\B6
			{
				if (iter->size() == k)//ѡ\D4\F1\C6\E4\D6б\DF\CA\FD\B5\C8\D3\DA
				{
					for (int j = 0; j < root.neighborIdx.size(); j++)
					{
						for (int index = 0; index < iter->size(); index++)
						{

							unsigned int s2 = 0;
							for (; s2 < root.neighborIdx.size(); s2++)
							{
								if (iter->at(index) == s2)
									break;

							}
							if (j > s2 && index == iter->size() - 1)
							{
								cutedpoint = *iter;
								//cutedpoint.push_back(root.neighborIdx.at(j));
								cutedpoint.push_back(j);//store the index of neighbor array
								root.cutedpoints.push_back(cutedpoint);
							}
						}
					}
				}
			}
		}
	}
	return root.cutedpoints.size();
}

bool Mesh3D::computeSubgraph(int flag)
{
	for (int j = 0; j < pvertices_list_->size(); j++)
	{
		if (pvertices_list_->at(j)->subgraphflag != flag)
		{
			continue;
		}
		else
		{
			// bian you  liangge dian 

			for (int k = 0; k < pvertices_list_->size(); k++)
			{
				pvertices_list_->at(k)->sumVector = (0.0, 0.0, 0.0);
				pvertices_list_->at(k)->set_seleted(UNSELECTED);
			}
			ComputeSumVector(j, flag);

			for (int k = 0; k < pvertices_list_->size(); k++)
			{
				//pvertices_list_.at(k).sumVector = (0.0, 0.0, 0.0);
				pvertices_list_->at(k)->truth =true;
				pvertices_list_->at(k)->set_seleted(UNSELECTED);
			}
			ComputeAngle(j, pvertices_list_->at(j)->sumVector);
			if (pvertices_list_->at(j)->truth)
			{
				return true;
			}
			for (int k = 0; k < pvertices_list_->size(); k++)
			{
				pvertices_list_->at(k)->sumVector = (0.0, 0.0, 0.0);
				pvertices_list_->at(k)->truth = true;
				pvertices_list_->at(k)->set_seleted(UNSELECTED);
			}
		}
	}
	return false;
}

void Mesh3D::ComputeSumVector(int id, int flag)
{
	HE_vert& root = *(pvertices_list_->at(id));
	root.set_seleted(SELECTED);
	for (auto iter = pvertices_list_->at(id)->neighborIdx.begin(); iter != pvertices_list_->at(id)->neighborIdx.end(); iter++)
	{
		HE_vert& neighbor_ =*(pvertices_list_->at(*iter));
		if (neighbor_.selected() == SELECTED)
		{
			continue;
		}
		else
		{
			Vec3f vector_ = neighbor_.position() - root.position();
			ComputeSumVector(*iter, flag);
			root.sumVector = root.sumVector + vector_ + neighbor_.sumVector;
		}
	}
	return;
}

void Mesh3D::ComputeAngle(int id, Vec3f sum)
{
	HE_vert& root = *(pvertices_list_->at(id));
	root.set_seleted(SELECTED);

	for (auto iter = pvertices_list_->at(id)->neighborIdx.begin(); iter != pvertices_list_->at(id)->neighborIdx.end(); iter++)
	{
		HE_vert& neighbor_ =* (pvertices_list_->at(*iter));
		if ((neighbor_.selected() == SELECTED))
		{
			continue;
		}
		else
		{
			ComputeAngle(neighbor_.id(), sum);
			Vec3f vector_ = neighbor_.position() - root.position();
			vector_.normalize();
			bool temp;
			if (vector_.dot(sum / sum.length()) > 0.3420)
			{
				temp = true;
			}
			else
			{
				temp = false;
			}
			root.truth = temp&& root.truth&&neighbor_.truth;
		}
	}
	return;
}

bool Mesh3D::findsubgraph(int * roots_, int num_roots_, int * cuts, int num_cuts_)
{
	
	for (int i = 0; i < num_roots_; i++)
	{
		//cout << "root: " << roots_[i] + 1;
		HE_vert& tosplitVert_ = *(pvertices_list_->at(roots_[i]));
		InsertVertex(tosplitVert_.position());
		(*(pvertices_list_->rbegin()))->vert_pair_ = roots_[i];
		for (auto iter = tosplitVert_.cutedpoints[cuts[i]].begin(); iter != tosplitVert_.cutedpoints[cuts[i]].end(); iter++)
		{
			(*(pvertices_list_->rbegin()))->neighbor_search_.push_back(tosplitVert_.neighborIdx[*iter]);//add neighbor to new vertex
			HE_vert* vert_neighbor_ = pvertices_list_->at(tosplitVert_.neighborIdx[*iter]);
			auto find_ = find(vert_neighbor_->neighbor_search_.begin(), vert_neighbor_->neighbor_search_.end(), tosplitVert_.id());//modify neighbor vertex's neighbor
			if (find_ != vert_neighbor_->neighbor_search_.end())
			{
				*find_ = pvertices_list_->size()-1;
				tosplitVert_.neighbor_search_.erase(find(tosplitVert_.neighbor_search_.begin(), tosplitVert_.neighbor_search_.end(), vert_neighbor_->id()));
			}
		}
		if (tosplitVert_.neighbor_search_.size()==0)
		{
			return false;
		}
	}
	for (int i = 0; i < pvertices_list_->size(); i++)
	{
		pvertices_list_->at(i)->set_seleted(UNSELECTED);
		pvertices_list_->at(i)->subgraphflag=-1;
	}
	int flag = -1;
 	for (int j = 0; j < pvertices_list_->size(); j++)
	{
		if (pvertices_list_->at(j)->selected() == UNSELECTED)
		{
			flag++;
			DFS((pvertices_list_->at(j)), (pvertices_list_->at(j)), flag);
		}
	}
	return true;
}

void Mesh3D::DFS(HE_vert* last_, HE_vert* cur, int flag)
{
	if (cur->selected() == UNSELECTED)
	{
		cur->subgraphflag = flag;
		cur->set_seleted(SELECTED);
		//cur->subgraphflag.push_back(flag);
		if (last_ != cur)
		{
			HE_vert* start_=last_;
			HE_vert* end_=cur;
			if (last_->id()>=num_vertex_)
			{
				start_ = pvertices_list_->at(last_->vert_pair_);
			}
			if (cur->id() >= num_vertex_)
			{
				end_ = pvertices_list_->at(cur->vert_pair_);
			}
			InsertEdge(start_,end_)->subgraphflag=flag;
		}
		for (int i = 0; i < cur->neighbor_search_.size(); i++)
		{
			last_ = cur;
			DFS(last_,pvertices_list_->at(cur->neighbor_search_[i]), flag);
		}
	}
}

void Mesh3D::Reset(int roots[],int num)
{
	for (int i = 0; i < num; i++)
	{
		pvertices_list_->at(roots[i])->neighbor_search_ = pvertices_list_->at(roots[i])->neighborIdx;
		pvertices_list_->pop_back();
	}
	//restore 
	for (int j = 0; j < pvertices_list_->size(); j++)
	{
		pvertices_list_->at(j)->set_seleted(UNSELECTED);
		pvertices_list_->at(j)->neighbor_search_ = pvertices_list_->at(j)->neighborIdx;
	}
}
