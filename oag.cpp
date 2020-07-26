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
typedef Arrangement_2::Face_handle                                      Face_handle;
typedef CGAL::Triangular_expansion_visibility_2<Arrangement_2>          TEV;

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

double* resolution(Polygon& p)
{
  double* dxy = new double[2];
  dxy[0] = 999999;
  dxy[1] = 999999;
  unsigned int n = p.size();
  double d;

  if (abs(p[n-1][0]-p[0][0]) > 0)
    dxy[0] = abs(p[n-1][0]-p[0][0]);

  if (abs(p[n-1][1]-p[0][1]) > 0)
    dxy[1] = abs(p[n-1][1]-p[0][1]);
  
  for (size_t i = 0; i < n-1; i++)
  {
    d = abs(p[i][0]-p[i+1][0]);
    if (d != 0 and d<dxy[0])
      dxy[0] = d;

    d = abs(p[i][1]-p[i+1][1]);
    if (d != 0 and d<dxy[1])
      dxy[1] = d;
  }
  
  return dxy;
}

vector<Point> discrate(Polygon& p)
{
  vector<Point> points;
  for (size_t i = 0; i < p.size(); i++)
    points.push_back(p[i]); 

  CGAL::Bbox_2 box = p.bbox();
  double x = box.xmin();
  double y = box.ymin();
  double xmax = box.xmax();
  double ymax = box.ymax();
  double* resolutionxy = resolution(p);

  while (x <= xmax)
  {
    y = box.ymin();
    while (y <= ymax)
    {
      
      if (p.is_simple())
      {
        Point point(x,y);
        int bs = p.bounded_side(point);

        if(bs==CGAL::ON_BOUNDED_SIDE || bs==CGAL::ON_BOUNDARY) //  || bs==CGAL::ON_BOUNDARY
          points.push_back(point);
      }
      y += resolutionxy[1];
    }
    x += resolutionxy[0];
  } 
  
  return points;
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

    vp.push_back(arrangement_to_polygon(regular_output));
  }

  return vp;
}

int** matrix_visibility(vector<Polygon> vp, int n, vector<Point> discrate, int m)
{

  int **arr = (int **)malloc(m * sizeof(int *));
  for (size_t i=0; i<m; i++)
    arr[i] = (int *)malloc(n * sizeof(int));
  
  for (size_t i = 0; i < n; i++)
  {
    //Polygon pi = visibility(points, n, i);
    for (size_t j = 0; j < m; j++)
    {
        int bs = vp[i].bounded_side(discrate[j]);
        
        if(bs==CGAL::ON_BOUNDED_SIDE || bs==CGAL::ON_BOUNDARY)
          arr[j][i] = 1;
          else
            arr[j][i] = 0;
          
    }
  }

  return arr;
}

int** matrix_visibility_refresh(int** matrix, vector<Polygon> vp, int n, vector<Point> dis, int m)
{
  int new_m = dis.size();
  int **arr = (int **)malloc(new_m * sizeof(int *));
  for (size_t i=0; i<new_m; i++)
    arr[i] = (int *)malloc(n * sizeof(int));

  for (size_t i = 0; i < m; i++)
    for (size_t j = 0; j < n; j++)
      arr[i][j] = matrix[i][j];
    
  for (size_t i = 0; i < n; i++)
    for (size_t j = m; j < new_m; j++)
    {
        int bs = vp[i].bounded_side(dis[j]);
        
        if(bs==CGAL::ON_BOUNDED_SIDE || bs==CGAL::ON_BOUNDARY)
          arr[j][i] = 1;
          else
            arr[j][i] = 0;
          
    }

  return arr;
}

void write_problem(string& filename, int** matrix, int& n, int& m, int& b)
{
  ofstream myfile (filename);
  if (myfile.is_open())
  {
    string cnf = "";
    int ind = 0;
    int num_var = n + n*(n-b);
    int num_cla =  m + n*(n-b) + (n-b) + (n*(n-b-1)*(n-b))/2;
    myfile<<"p cnf " + to_string(num_var) + " " + to_string(num_cla) + "\n";

    for (size_t i = 0; i < n; i++)
      for (size_t k = 0; k < n-b; k++)
      {
        ind = -b*i;
        cnf = "-" + to_string((n + 1 + (i*n+k))+ind) + " -" + to_string((i+1)) + " 0\n";
        myfile<<cnf;
      }

    for (size_t k = 0; k < n-b; k++)
    {
      cnf = "";
      for (size_t i = 0; i < n; i++)
      {
        ind = -b*i;
        cnf += to_string((n + 1 + (i*n+k))+ind) + " ";
      }
      cnf += "0\n";
      myfile<<cnf;
    }

    for (size_t i = 0; i < n; i++)
      for (size_t t = 0; t < n-b; t++)
        for (size_t k = 0; k < t; k++)
          {
            ind = -b*i;
            cnf = "-" + to_string((n + 1 + (i*n+k))+ind) + " -" + to_string((n + 1 + (i*n+t))+ind) + " 0\n";
            myfile<<cnf;
          }

    for (size_t i = 0; i < m; i++)
    {
      cnf = "";
      for (size_t j = 0; j < n; j++)
        if(matrix[i][j]!=0)
          cnf += to_string(j+1) + " ";
      
      cnf += "0\n";
      myfile<<cnf;
    }
    
    myfile.close();
  }
  else cout << "Unable to open file";
}

vector<int> read_solve(string& filename, int& n)
{
  vector<int> result;

  string text;
  ifstream MyReadFile(filename);
  int status = 0, i = 0;

  while(getline (MyReadFile, text))
  {
    if(i == 0)
      if(text.compare("SAT")==0)
        status = 1;
    
    if(i==1 && status == 1)
    {
      vector<string> line = split(text, " ");
      for (size_t j = 0; j < n; j++)
        if(line[j].find("-") == std::string::npos)
          result.push_back(std::stoi(line[j]) - 1);
    }
    i++;
  }
  MyReadFile.close(); 


  return result;
}

vector<int> solver(string& filename, int** a, int& n, int& m)
{
  vector<int> result;
  string cnf_out = "result/cnf/" + filename + ".txt";
  string file_solve = "result/cnf/solve/" + filename + ".txt";
  string time_limit = std::to_string(n/8);
  string str = "minisat -cpu-lim=" + time_limit + " " + cnf_out + " " + file_solve;
  const char *command = str.c_str();
  int b = n/4;

  while(1)
  {
    cout<<"==========================================================="<<endl;
    cout<<"                   BROJ CUVARA  >   "<< b <<"     <"<<endl;
    cout<<"==========================================================="<<endl;

    // write problem
    write_problem(cnf_out, a, n, m, b);

    // solving problem
    system(command);

    // read solve 
    vector<int> onps = read_solve(file_solve, n);
    
    // if onps is emptey then SOLVE is INDETERMINATE or UNSAT
    b = onps.size()-1;
    if(b == -1) break;
    else result = onps;
  }
  
  return result;
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

  outfile.open("test03-tred1.txt", std::ios_base::app); // append instead of overwrite
  outfile << tekst; 
}


int main(int argc, char const *argv[])
{
  for(int file_number=198; file_number<199; file_number+=2) // start = 8
  for(int file_order=72; file_order<=file_number; file_order++) // file_order<=file_number
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

    // discretisation
    cout<<"Discretisation............................";
    auto start = high_resolution_clock::now(); 
    vector<Point> d = discrate(p);
    int m = d.size();
    cout<<"end!!!"<<endl;

    // visibility polygons
    cout<<"Creating visibility polygons..............";
    vector<Polygon> vp = visibility_polygon(points, n);
    cout<<"end!!!"<<endl;

    // create matrix visability
    int **a = matrix_visibility(vp, n, d, m);
    cout<<"Created matrix visability"<<endl;

    // pretrocesing time
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start); 
    cout << "Time Execution Pretrocesing >: " << duration.count() << " microseconds" << endl;
    write_test(std::to_string(duration.count()) + ";"); 

    int broj_ponavljanja = 0;
    int broj_dodanih_tacaka = 0;
    int broj_cuvara = n/4;
    start = high_resolution_clock::now(); 
    while(1)
    {
      // solving
      vector<int> s = solver(filename, a, n, m);
  
      // convert polygons - K to Kernel
      Point_2* points_2 = new Point_2[n];
      for (size_t i = 0; i < n; i++)
        points_2[i] = converter2(points[i]);//Point_2(points[i][0], points[i][1]);
      Polygon_2 p_2 = Polygon_2(points_2, points_2 + n);
      vector<Polygon_2> vp_2 = convertVP(vp);
      // converted

      //finding the region
      Polygon_set_2 S(p_2);
      for(size_t i = 0; i < s.size(); i++)
      {
        if(vp_2[s[i]].orientation()==-1)
          vp_2[s[i]].reverse_orientation();
        S.difference(vp_2[s[i]]);//cout<<points[s[i]]<<", ";
      }
      std::list<Polygon_with_holes_2> res;
      S.polygons_with_holes (std::back_inserter (res));

      // if there is no region then the end
      if(res.size()==0){
        write_test(std::to_string(s.size()) + ";" + std::to_string(broj_ponavljanja) + ";" + std::to_string(broj_dodanih_tacaka) + ";"); 
        broj_cuvara = s.size();
        break;
      } 

      // finding the centroid
      Pwh_list_2::const_iterator it_res;
      for (it_res = res.begin(); it_res != res.end(); ++it_res) {
        std::cout << "-->";
        // print_polygon_with_holes (*it_res);
        vector<Point_2> t = polygon_with_holesTOvector(*it_res);
        Point_2 c2 = CGAL::centroid(t.begin(), t.end(),CGAL::Dimension_tag<0>());
        d.push_back(converter(c2));
        cout << c2 << endl;
        broj_dodanih_tacaka++;
      }
      if(d.size()>1 && d[d.size()-1]==d[d.size()-2]){
        write_test(std::to_string(s.size()) + ";" + std::to_string(broj_ponavljanja-1) + ";" + std::to_string(broj_dodanih_tacaka-1) + ";"); 
        break;
      }
      // reset problem
      for (size_t i = 0; i < s.size(); i++)
        cout<<s[i]<<", ";
      cout<<endl;
      
      //a = matrix_visibility(vp, n, d, m);
      int** b = matrix_visibility_refresh(a, vp, n, d, m);
      free(a);
      a = b;
      m = d.size();
      broj_ponavljanja++;
    }
    stop = high_resolution_clock::now();
    duration = duration_cast<microseconds>(stop - start); 
    cout << "Time Execution Pretrocesing >: " << duration.count() << " microseconds" << endl;
    write_test(std::to_string(duration.count()) + "\n");
    free(a);
  }

  return EXIT_SUCCESS;
}