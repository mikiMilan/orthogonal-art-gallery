#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_polygon_visibility_2.h>
#include <CGAL/Triangular_expansion_visibility_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Gps_segment_traits_2.h>
#include <CGAL/centroid.h>
#include <list>
#include <minisat/core/Solver.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <chrono>
#include <bits/stdc++.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel             K;
typedef CGAL::Polygon_2<K>                                              Polygon;
typedef CGAL::Point_2<K>                                                Point;
typedef CGAL::Exact_predicates_exact_constructions_kernel               Kernel;
typedef Kernel::Point_2                                                 Point_2;
typedef CGAL::Polygon_2<Kernel>                                         Polygon_2;
typedef Kernel::Segment_2                                               Segment_2;
typedef CGAL::Arr_segment_traits_2<Kernel>                              Traits_2;
typedef CGAL::Arrangement_2<Traits_2>                                   Arrangement_2;
typedef Arrangement_2::Edge_const_iterator                              Edge_const_iterator;
typedef Arrangement_2::Halfedge_const_handle                            Halfedge_const_handle;
typedef Arrangement_2::Halfedge_handle                                  Halfedge_handle;
typedef Arrangement_2::Halfedge_handle                                  Halfedge_handle;
typedef CGAL::Simple_polygon_visibility_2<Arrangement_2,CGAL::Tag_true> RSPV;
typedef CGAL::Polygon_with_holes_2<Kernel>                              Polygon_with_holes_2;
typedef std::list<Polygon_with_holes_2>                                 Pwh_list_2;
typedef CGAL::Gps_segment_traits_2<Kernel>                              Traits_22;
typedef CGAL::General_polygon_set_2<Traits_22>                          Polygon_set_2;

typedef Arrangement_2::Face_handle                              Face_handle;
// Define the used visibility class 
typedef CGAL::Triangular_expansion_visibility_2<Arrangement_2>  TEV;

CGAL::Cartesian_converter<Kernel,K> converter; 
CGAL::Cartesian_converter<K, Kernel> converter2; 

using namespace std;
using namespace std::chrono;

vector<string> split(const string& str, const string& delim)
{
    vector<string> tokens;
    size_t prev = 0, pos = 0;
    do
    {
        pos = str.find(delim, prev);
        if (pos == string::npos) pos = str.length();
        string token = str.substr(prev, pos-prev);
        if (!token.empty()) tokens.push_back(token);
        prev = pos + delim.length();
    }
    while (pos < str.length() && prev < str.length());

    return tokens;
}

Point* read_point_of_polygon(string& filename, int *n)
{
  Point* points;
  *n = 0;

  string text;
  ifstream MyReadFile(filename);

  while(getline (MyReadFile, text))
  {
    vector<string> line = split(text, " ");
    stringstream geek(line[0]); 
    *n = 0; 
    geek >> *n;
    points = new Point[*n];
    
    int numerator = 0, denominator = 1;
    vector<string> number;
    for (int i = 1; i < line.size(); i+=2)
    {
      number = split(line[i], "/");
      numerator = 0;
      denominator = 1;

      stringstream geekN(number[0]); 
      geekN >> numerator;
      stringstream geekD(number[1]);
      geekD >> denominator;

      double po1 =(double)numerator/denominator;

      number = split(line[i+1], "/");
      numerator = 0;
      denominator = 1;

      stringstream geekN2(number[0]); 
      geekN2 >> numerator;
      stringstream geekD2(number[1]);
      geekD2 >> denominator;

      double po2 =(double)numerator/(double)denominator;

      points[i/2] = Point(po1, po2);
    }

  }
  MyReadFile.close(); 

  return points;
}

Polygon read_polygon(string filename, int *n)
{
  Point* points = read_point_of_polygon(filename, n);
  return Polygon(points, points+*n);
}

Polygon arrangement_to_polygon(Arrangement_2& arr)
{
  Edge_const_iterator eit = arr.edges_begin();
  vector<Point_2> points_2 = {eit->source()->point(), eit->target()->point()};

  vector<Edge_const_iterator> edges_it;
  for(++eit; eit != arr.edges_end(); ++eit)
    edges_it.push_back(eit);

  while(edges_it.size()>0)
  {
    for (size_t i = 0; i < edges_it.size(); i++)
    {
      eit = edges_it[i];

      if(points_2[0]==eit->source()->point())
      {
        points_2.insert(points_2.begin(), eit->target()->point());
        edges_it.erase (edges_it.begin()+i);
        i--;
        continue;
      }
      if(points_2[0]==eit->target()->point())
      {
        points_2.insert(points_2.begin(), eit->source()->point());
        edges_it.erase (edges_it.begin()+i);
        i--;
        continue;
      }
      if(points_2[points_2.size()-1]==eit->source()->point())
      {
        points_2.push_back(eit->target()->point());
        edges_it.erase (edges_it.begin()+i);
        i--;
        continue;
      }
      if(points_2[points_2.size()-1]==eit->target()->point())
      {
        points_2.push_back(eit->source()->point());
        edges_it.erase (edges_it.begin()+i);
        i--;
      }
    }
  }
  

  int n = arr.number_of_edges();
  Point* points = new Point[n];

  for (size_t i = 0; i < n; i++)
    points[i] = converter(points_2[i]);
  
  return Polygon(points, points + n);
}

vector<Polygon> visibility_polygon(Point* p, int& n)
{
  vector<Polygon> vp;
  
  // convert to exact computation paradigm
  Point_2* points = new Point_2[n];
  for (size_t i = 0; i < n; i++)
    points[i] = Point_2(p[i][0], p[i][1]);
  // end - convert to exact computation paradigm

  std::vector<Segment_2> segments;
  for (size_t i = 0; i < n-1; i++)
    segments.push_back(Segment_2(points[i], points[i+1]));
  segments.push_back(Segment_2(points[n-1], points[0]));

  Arrangement_2 env;
  CGAL::insert_non_intersecting_curves(env,segments.begin(),segments.end());
  Arrangement_2 regular_output;

  for (size_t i = 0; i < n; i++)
  {
    
    Halfedge_const_handle he = env.halfedges_begin();
    int q_1 = i-1;
    if (q_1==-1)
      q_1 = n-1;
    
    while (he->source()->point() != points[q_1] || he->target()->point() != points[i])
      he++;

    TEV tev(env);
    Face_handle fh = tev.compute_visibility(points[i], he, regular_output);
    /*
    RSPV regular_visibility(env);
    regular_visibility.compute_visibility(points[i], he, regular_output);
    */

    vp.push_back(arrangement_to_polygon(regular_output));
  }

  return vp;
}

template<class Kernel, class Container>
void print_polygon (const CGAL::Polygon_2<Kernel, Container>& P)
{
  typename CGAL::Polygon_2<Kernel, Container>::Vertex_const_iterator vit;
  std::cout << "[ " << P.size() << " vertices:";
  for (vit = P.vertices_begin(); vit != P.vertices_end(); ++vit)
    std::cout << " (" << *vit << ')';
  std::cout << " ]" << std::endl;
}

template<class Kernel, class Container>
void print_polygon_with_holes(const CGAL::Polygon_with_holes_2<Kernel, Container> & pwh)
{
  if (! pwh.is_unbounded()) {
    std::cout << "{ Outer boundary = "; 
    print_polygon (pwh.outer_boundary());
  }
}

void print_set_polygons(Polygon_set_2 S)
{
  std::list<Polygon_with_holes_2> res;
  S.polygons_with_holes (std::back_inserter (res));

  Pwh_list_2::const_iterator it_res;
  for (it_res = res.begin(); it_res != res.end(); ++it_res) 
    print_polygon_with_holes(*it_res);
}

template<class Kernel, class Container>
Traits_2::FT area_polygon_with_holes(const CGAL::Polygon_with_holes_2<Kernel, Container> & pwh)
{
  Traits_2::FT res = 0;
  if (! pwh.is_unbounded()) {
    res = pwh.outer_boundary().area();
  } else
    return res;

  typename CGAL::Polygon_with_holes_2<Kernel,Container>::Hole_const_iterator hit;
  unsigned int k = 1;
  for (hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit, ++k) 
    res -= (*hit).area();

  return res;
}

Traits_2::FT area_set_polygons(Polygon_set_2 S)
{
  Traits_2::FT povrsina = 0;

  std::list<Polygon_with_holes_2> res;
  S.polygons_with_holes (std::back_inserter (res));

  Pwh_list_2::const_iterator it_res;
  for (it_res = res.begin(); it_res != res.end(); ++it_res) 
    povrsina += area_polygon_with_holes(*it_res);
      
  return povrsina;
}

template<class Kernel, class Container>
vector<Point_2> polygon_with_holesTOvector(const CGAL::Polygon_with_holes_2<Kernel, Container> & pwh)
{
  vector<Point_2> pol;

  if (! pwh.is_unbounded())
  {
    Polygon_2 P = pwh.outer_boundary();
    if(P.is_simple())
    {
      typename CGAL::Polygon_2<Kernel, Container>::Vertex_const_iterator vit;
      for (vit = P.vertices_begin(); vit != P.vertices_end(); ++vit)
        pol.push_back(*vit);
    }else{
      typename CGAL::Polygon_2<Kernel, Container>::Vertex_const_iterator i, j, start, end;
      bool uslov = 1;
      for (i = P.vertices_begin(); i != P.vertices_end() && uslov; ++i)
        for (j = i+1; j != P.vertices_end() && uslov; ++j)
          if(*i==*j)
          {
            start = i;
            end = j;
            uslov = 0;
          }
      for(i = start; i != end; ++i)
        pol.push_back(*i);
    }
  }

  return pol;
}

vector<Polygon_2> convertVP(vector<Polygon> vp)
{
  vector<Polygon_2> vp_2;

  for (size_t i = 0; i < vp.size(); i++)
  {
    Point_2* p2 = new Point_2[vp[i].size()];

    for (size_t j = 0; j < vp[i].size(); j++)
      p2[j] = converter2(vp[i][j]);
    
    vp_2.push_back(Polygon_2(p2, p2 + vp[i].size()));
  }
  
  return vp_2;
}

void write_test(string tekst)
{
  std::ofstream outfile;

  outfile.open("test-rand-4.txt", std::ios_base::app); // append instead of overwrite
  outfile << tekst; 
}

bool emptyforeach(vector<Polygon_set_2> vp)
{
  for (size_t i = 0; i < vp.size(); i++)
  {
    std::list<Polygon_with_holes_2> res;
    vp[i].polygons_with_holes(std::back_inserter(res));

    // if there is no region then the end
    if(res.size()!=0)
      return 0;
  }
  
  return 1;
}

int main(int argc, char const *argv[])
{
  for(int file_number=186; file_number<201; file_number+=2) // start = 8
  for(int file_order=1; file_order<=file_number; file_order++) // file_order<=file_number
  {
    string filename = "rand-" + std::to_string(file_number) + "-" + std::to_string(file_order) + ".pol";

    write_test(std::to_string(file_number) + "-" + std::to_string(file_order) + ";");
    //string filename = "rand-200-2.pol";
    string location = "instance/random/" + filename;

    // read polygon
    int n;
    Point* points = read_point_of_polygon(location, &n);
    Polygon p = Polygon(points, points+n);
    cout<<"Created polygon "<<filename<<endl;

    auto start = high_resolution_clock::now();
    // visibility polygons
    cout<<"Creating visibility polygons..............";
    vector<Polygon> vp = visibility_polygon(points, n);
    //vector<Polygon_2> vp_2 = convertVP(vp); //mozda zatreba
    //points_2[i] = converter2(points[i]);//Point_2(points[i][0], points[i][1]);
    cout<<"end!!!"<<endl;

    // convert polygons - K to Kernel
    Point_2* points_2 = new Point_2[n];
    for (size_t i = 0; i < n; i++)
      points_2[i] = converter2(points[i]);//Point_2(points[i][0], points[i][1]);
    Polygon_2 p_2 = Polygon_2(points_2, points_2 + n);
      
    // convert visability to vp_2 and vp_S
    vector<Polygon_2> vp_2 = convertVP(vp);

    for(size_t i = 0; i < vp_2.size(); i++)
      if(vp_2[i].orientation()==-1)
        vp_2[i].reverse_orientation();
    
    vector<Polygon_set_2> vp_S;
    for (size_t i = 0; i < vp_2.size(); i++)
      vp_S.push_back(Polygon_set_2(vp_2[i]));
    // converted

    vector<int> I;

    cout<<"Trazim cuvare...."<<endl;
    while(!emptyforeach(vp_S))
    {
      Polygon_set_2 max_vp = vp_S[0];
      int max_indeks = 0;

      for (size_t i = 1; i < vp_S.size(); i++)
        if(area_set_polygons(vp_S[i]) > area_set_polygons(max_vp))
        {
          max_vp = vp_S[i];
          max_indeks = i;
        }
        
      I.push_back(max_indeks);
      cout<<p[max_indeks]<<endl;

      // smanjenje povrsine ostalim regijama vp
      for (size_t i = 0; i < vp_S.size(); i++)
      vp_S[i].difference(max_vp);
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start); 
    cout << "Time Execution>: " << duration.count() << " microseconds" << endl;
    write_test(std::to_string(I.size()) + ";" + std::to_string(duration.count()) + "\n");
    cout<<"Dovoljan broj cuvara je: "<<I.size()<<endl;

  }

  return EXIT_SUCCESS;
}