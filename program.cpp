#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <string.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <sys/poll.h>
#include <netinet/in.h>

using namespace std;

#define PI 3.14159265
#define eps 0.0000001


ofstream LOG("data/log.txt", std::ios::app);

void writeToClient(int clientSocket, stringstream& strOut) {
	string msg = strOut.str();
	LOG << "TO [" << clientSocket << "]: " << msg << endl;
	send(clientSocket, msg.c_str(), msg.size() + 1, 0);
	strOut.str("");
}

struct  Node {
	double aX = -1; 
    double aY = -1;
    Node *l; 
    Node *r; 
};

class Point {
	// класс точка
	private:
		double x, y; // coordinates
	public:
		int cluster;
		Point(double xx, double yy) : x(xx), y(yy), cluster(-1) {};
        ~Point() {}
       // Point(const Point& p) : x(p.x), y(p.y), cluster(p.cluster) {}
		double getX() {
			return x;
		}
		double getY() {
			return y;
		}
		double distanceTo(Point p) {
		// расстояние от этой точки до точки (xx, yy)
			return sqrt((p.getX() - x)*(p.getX()-x) + (p.getY()-y)*(p.getY()-y));
		}
		Point& operator =(const Point& p)  {x=p.x;y=p.y; cluster=p.cluster; return *this;}
		void operator +=(const Point& p) {x+=p.x;y+=p.y;}
		void operator *=(const Point& p) {x*=p.x;y*=p.y;}
		Node *Tree = NULL; // для дендрограммы
};


class Matrix2D {
	private:
		vector<vector<double> > matrix;
	public:
		Matrix2D() {
			matrix.resize(2);
    		for (unsigned i = 0; i < matrix.size(); i++) {
        		matrix[i].resize(2);
    		}
		};
		void set(int i, int j, double x) {
			matrix[i][j] = x;
		}
		double get(int i, int j) {
			return matrix[i][j];
		}
		void make_covariance(vector<Point>& points, Point center) {
			// Calculating covariance matrix of points
			double value, temp1, temp2;
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					value = 0;
					for (size_t k = 0; k < points.size(); k++) {
						if (i == 0) temp1 = points[k].getX() - center.getX();
						else temp1 = points[k].getY() - center.getY();

						if (j == 0) temp2 = points[k].getX() - center.getX();
						else temp2 = points[k].getY() - center.getY();
						value+=temp1*temp2;
					}
					value/=points.size();
					matrix[i][j] = value;
				}
			}

		}
		vector <double> eighensystem() {
			// возвращаем массив собственных векторов и значений
			vector <double> values(2); // собственные значения
			vector <double> vectors(6); // 2 вектора и 2 значения 
			double D; // Discriminant
			D = (matrix[0][0]+matrix[1][1])*(matrix[0][0]+matrix[1][1]) - 
					4*(matrix[0][0]*matrix[1][1]-matrix[0][1]*matrix[1][0]);
			//adding roots
			values[0] = (matrix[0][0]+matrix[1][1]+ sqrt(D))/2;
			values[1] = (matrix[0][0]+matrix[1][1]- sqrt(D))/2;
			for (int i = 0; i < 2; i++) {
				vectors[2*i] = matrix[0][1]/(values[i] - matrix[0][0]); // vectors[i+1] = 1;
                vectors[2*i+1] = 2*sqrt(values[i])/sqrt(vectors[2*i]*vectors[2*i] + 1);
				vectors[2*i] *= vectors[2*i+1]; 
			}
			vectors[4] = values[0];
			vectors[5] = values[1];
			return vectors;
		}
};

class Group {
	private:
		// массив точек, содержащихся в группе
		vector <Point> list;
		// число точек в группе
		int p_count; 
		// центра масс
		Point gr_center;
		// цвет группы
		int r = rand() % 256, g = rand() % 256, b = rand() % 256;
	public:
		Group (double x, double y, double disp1, double disp2, int N) : list(), p_count(N), gr_center(0,0){
			// Генерация N - точек
			for (int i = 0; i < N; i++) {
				double s1 = 0, s2 = 0;
				for (int j = 0; j < 1000; j++) {
					s1+= (rand() % 1000) * 0.002 - 1;
					s2+= (rand() % 1000) * 0.002 - 1;
				}
				Point p(x + disp1* s1/1000, y + disp2 * s2/1000);	
				list.push_back(p);
			}
			p_count = list.size();
			if (p_count != 0) gr_center = mCenter();

		}
		void add_point(Point p) { 
		// добавление точки в группу
			list.push_back(p);
			p_count++;
			gr_center = mCenter();
		}
		int get_size() {
			// вывод количества точек
			return list.size();
		}

		Point mCenter() {
			// Возвращаем центр масс группы
			double xs = 0, ys = 0;
			for (int i = 0; i < p_count; i++) {
				xs+= list[i].getX();
				ys+= list[i].getY();
			}
			Point p(xs/p_count, ys/p_count);
			return p;
		}
		void move(double x, double y) {
			// смещение группы на вектор {x, y}
			Point p(x, y);
			for (int i = 0; i < p_count; i++)
				list[i]+= p;
			gr_center = mCenter();
		}

		void rotate(double alfa) {
			// Вращения отн-но центра масс на alfa по часовой
			double x, y;
			Point p = mCenter();
			for (int i = 0; i < p_count; i++) {
				x = list[i].getX() - p.getX();
				y = list[i].getY() - p.getY();
				list[i] = Point(x*cos(alfa*PI/180) + y*sin(alfa*PI/180) + p.getX(),
								y*cos(alfa*PI/180) - x*sin(alfa*PI/180) + p.getY());
			}
		}
		void scale (double k1, double k2) {
			// растяжение группы в k1 - раз по Х, k2 - раз по Y
			Point p(k1, k2);
			for (int i = 0; i < p_count; i++)
				list[i]*= p;
			gr_center = mCenter();
		}
		void print_group(ofstream& print, bool colourful) { 
		// вывод всех точек группы
			for (int i = 0; i < p_count; i++) {
				print << list[i].getX() << " " << list[i].getY();
				if (colourful) {
					print << " " << r << " " << g << " " << b;
				}
				print << endl;
			}
		}
		friend class Field;
};

class Cluster {
	private:
		vector <Point> points;
		int size = 0;
		int r = rand() % 256, g = rand() % 256, b = rand() % 256;
	public:
		Point c{0, 0};
		Cluster(Point p) : points(), c(p) {}
		void add_point(Point p) {
			points.push_back(p);
			size++;
		}
		Point mCenter() {
			// returns center of a cluster
			double xs = 0, ys = 0;
			for (int i = 0; i < size; i++) {
				xs+= points[i].getX();
				ys+= points[i].getY();
			}
			Point p(1,1);
			if (size != 0) {
                Point t(xs/size, ys/size);
                p *= t;
            }
			return p;
		}
		void print_arrow(ofstream& print) {
			// отрисовка компонент кластера
			Matrix2D cov; // ковариоционная матрица
			Point newCenter = mCenter(); 
			vector <double> eighen; // Собственные вектора матрицы
			double value = 0;
			double temp1 = 0, temp2 = 0;
            if (points.size() < 10) return;
			// Заполняем ковариоцинную матрицу
			cov.make_covariance(points, newCenter);

			eighen = cov.eighensystem();
			// Печатаем
			print << newCenter.getX() << " " << newCenter.getY() << " " << 
			eighen[0] << " " << eighen[1] << endl;
			print << newCenter.getX() << " " << newCenter.getY() << " " << 
			eighen[2] << " " << eighen[3] << endl; 

		}
		void print_cluster(ofstream& print) { 
		// вывод всех точек группы
			for (int i = 0; i < (int)points.size(); i++) 
				print << points[i].getX() << " " << points[i].getY() << " " << r << " " << g << " " << b << endl;
	
		}
		void clear() {
			// удалить все точки кластера
			size = 0;
			points.clear();
		}
		friend class Field;
		friend class EM_alghoritm;
};

class Cluster_search {
	public:
		vector <Cluster> cl_list; // найденные кластеры
		vector <string> params; // параметры поиска
		Cluster_search() : cl_list() {}

		Cluster& get_cluster(int cluster_id) {
			return cl_list[cluster_id];
		}
		int clusterCount() {
			return cl_list.size();
		}
		string getParameter(int i) {
			if (i >= 0 && i < params.size())
				return params[i];
			else 
				return "";
		}
		void deleteClusters() {
			for (size_t i = 0; i < cl_list.size(); i++)
				cl_list[i].clear();
		}
		void print_pca() {
			// печать факторов
			ofstream print("data/pca.txt");
			for (int i = 0; i < (int)cl_list.size(); i++) {
				cl_list[i].print_arrow(print);
			}
		}
		void print_all(ofstream& print) {
			// печать всех кластеров
			for (int i = 0; i < (int)cl_list.size(); i++) {
				cl_list[i].print_cluster(print);
			}
		}
		friend class Field;
};

class Field {
	// класс поле
	private:
		int size = 0; // количество точек
		int sCount = -1; // номер последнего поиска
		bool ClusterMode = false; // можно ли запускать поиски кластеров 
		vector <Point *> points; // указатели на точки групп
		vector <Cluster_search> searches; 
	public:
		vector <Group> gr_list;
		double **distances; // матрица расстояний

		Point& get_point(int id) {
			// Возвращаем точку по ее id
			return *points[id];
		}
		bool getMode() {
			// true - we cluster, false - generating
			return ClusterMode;
		}
		void changeMode() {
			if (ClusterMode == false)
				ClusterMode = true;
			else
				ClusterMode = false;
		}
		int count() {
			// количество точек в поле
			return size;
		}
		void add_search(Cluster_search newSearch) {
			// добавить новый поиск кластеров
			searches.push_back(newSearch);
			sCount++;
		}
		int search_count() {
			// Получить  номер последнегоп оиска
			return sCount;
		}
		Cluster_search get_search(int search_id) {
			// Возвращаем поиск по его номеру
			if (search_id < 0 || search_id > sCount)
				return searches[sCount];
			else
				return searches[search_id];
		} 

		void add(Group g) {
			// добавляем группу в поле
			gr_list.push_back(g);
			for (int i = 0; i < g.get_size();i++) {
				points.push_back(new Point(2,3));
				points[size+i] = &gr_list[gr_list.size() - 1].list[i];
			}
			size+=g.p_count;
		}
		void clear() {
			// Очистить поле
			gr_list.clear();
			points.clear();
			searches.clear();
			size = 0;
			sCount = -1;
		}

		void load (ifstream& file) {
			double x, y;
			int id = 0, cnt = 0;
			string pr;
			int r1 = -1, g1 = -1, b1 = -1, rc, gc, bc;
			file >> cnt;
			// Загрузка групп
			for (int i = 0; i < cnt; i++) {
				file >> x >> y >> rc >> gc >> bc;
				if (rc != r1 || gc != g1 || bc != b1) {
					Group gr(0, 0, 0, 0, 0); // ! 
					add(gr);
					id = gr_list.size()-1;
				}
				Point p(x, y);
				gr_list[id].add_point(p);
				size++;
				r1 = rc;
				g1 = gc;
				b1 = bc;
			}
			// Загрузка кластеров
			while (file >> cnt) {
				getline(file, pr);
				getline(file, pr);
				r1 = g1 = b1 = -1;
				Cluster_search new_search;

				size_t pos = 0; // Считываем параметры
				while ((pos = pr.find(" ")) != string::npos) {
   				 	new_search.params.push_back(pr.substr(0, pos));
   					 pr.erase(0, pos + 1);
				}
				new_search.params.push_back(pr);

				for (int i = 0; i < cnt; i++) {
					file >> x >> y >> rc >> gc >> bc;
					Point p(x, y);
					if (rc != r1 || gc != g1 || bc != b1) {
						Cluster my_cluster(p);
						new_search.cl_list.push_back(my_cluster);
						id = new_search.cl_list.size()-1;
					}
					new_search.cl_list[id].add_point(p);
					r1 = rc;
					g1 = gc;
					b1 = bc;
				}
				searches.push_back(new_search);
				sCount++;
			}
			cout << "field: loaded | Total point count: " << size << " | Total search count: " << sCount + 1 << endl;
            LOG << "field: loaded | Total point count: " << size << " | Total search count: " << sCount + 1 << endl;

		}

		void archieve() {
			// We save groups and cluster searches to file
			ofstream a("data/ARCHIEVE.txt");
			a << size << endl; // кол-во точек
			for (int i = 0; i < (int)gr_list.size(); i++)
				gr_list[i].print_group(a, true);
			for (int j = 0; j <= sCount; j++) {
				a << size << endl;
				for (auto pr : searches[j].params) a << pr << " ";
				a << endl;
				for (int i = 0; i < (int)searches[j].cl_list.size(); i++) {
					searches[j].cl_list[i].print_cluster(a);
				}
			}
			cout << "archieved to data/ARCHIEVE.txt" << endl;
            LOG << "archieved to data/ARCHIEVE.txt" << endl;
		}

		void fill_matrix() {
			// заполняем матрицу расстояниями
			distances = new double*[size]; 
   			for (int i = 0; i < size; i++) distances[i] = new double[size];

   			for (int i = 0; i < size;i++) {
				for (int j = i + 1; j < size; j++) {
					distances[i][j] = distances[j][i] = points[i]->distanceTo(*points[j]); 
				}
			}
		}
		Point random_point() {
			// возвращаем случайную точку поля
			return get_point(rand() % size);
		}
		void group_list() {
			for (int i = 0; i < (int)gr_list.size(); i++) {
				cout << i << " " << gr_list[i].p_count << " points(" << gr_list[i].gr_center.getX() << ", " << gr_list[i].gr_center.getY() << ")" << endl;
                LOG << i << " " << gr_list[i].p_count << " points(" << gr_list[i].gr_center.getX() << ", " << gr_list[i].gr_center.getY() << ")" << endl;
			}
		}
		void search_list() {
			for (int i = 0; i <= sCount; i++) {
				cout << "[" << i << "] params = { ";
                LOG << "[" << i << "] params = { ";
				for (auto pr : searches[i].params) {
                    cout << pr << ", ";
                    LOG << pr << ", ";
                }
				cout << "... }" << endl;
                LOG << "... }" << endl;
			}
		}
		void print_all() {
			// печать всех групп поля
			ofstream print("data/field.txt");
			for (int i = 0; i < (int)gr_list.size(); i++) {
				gr_list[i].print_group(print, true);
			}
		}
};

class Porog {
	private:
		Field& my_field;
		Cluster_search my_search;
		char** matrix; // граф
	public:
		Porog(Field& f) : my_field(f) {};
		void dfs(int pp, int step, vector <bool>& visit) {
			// pp - посещаемая вершина, step - компонента связности (id кластера)
			visit[pp] = true;
			my_search.cl_list[step].add_point(my_field.get_point(pp));
			for (int j = 0; j < my_field.count(); j++) {
				if (visit[j] == false && matrix[pp][j] == 1) {
					dfs(j, step, visit);
				}
			}
		}
		int start(double value) {
			// value - порог
			my_search.params.push_back("porog");
			my_search.params.push_back(to_string(value));
			vector <bool> visit(my_field.count(), false);
			my_field.fill_matrix();
			// Матрица смежности графа
			matrix = new char*[my_field.count()]; 
   			for (int i = 0; i < my_field.count(); i++) 
        		matrix[i] = new char[my_field.count()];

			for (int i = 0; i < my_field.count();i++) { 
				for (int j = 0; j < my_field.count();j++) {
					if (my_field.distances[i][j] < value)
						matrix[i][j] = true;
					else 
						matrix[i][j] = false;
				}
			}

			int step = 0, pp = 0;

			while (pp < my_field.count()) { // пока не посетили все точки
				Cluster c(my_field.random_point());
				my_search.cl_list.push_back(c);
				dfs(pp, step, visit);
				for (int i = 0; i <= my_field.count(); i++) {
					if (visit[i] == false) {
						pp = i;
						break;
					}
				}
				step++;
			}
			cout << step << " clusters found" << endl;
			my_field.add_search(my_search);
			return my_field.search_count(); // return search_id
		}

		bool print(int search_id, ofstream& print) {
			my_search = my_field.get_search(search_id);
			my_search.print_all(print);
			FILE *f;
			f = popen("gnuplot kmeans.script", "r");
			return pclose(f);
		}
};

class Kmeans {
	private:
		Field& my_field;
		Cluster_search my_search;
	public:
		Kmeans(Field& f) : my_field(f) {};
		bool fill_clusters(int k) {
			// заполняем кластеры точками поля
			bool is_changed = false;
			for (int i = 0; i < my_field.count(); i++) {
				Point p = my_field.get_point(i);
				int candidate = 0;
				double distance = p.distanceTo(my_search.cl_list[0].c);
				for (int l = 1; l < k;l++) {
					if (p.distanceTo(my_search.cl_list[l].c) < distance) {
						distance = p.distanceTo(my_search.cl_list[l].c);
						candidate = l;
					}
				}
				if (p.cluster != candidate) {
					is_changed = true;
					my_field.get_point(i).cluster = candidate;
				}
				my_search.cl_list[candidate].add_point(p);
			}
			return is_changed;
		}

		int start(int k) { // алгоритм k-средних
			my_search.params.push_back("k-means");
			my_search.params.push_back(to_string(k));
			// Выбираем случайные k точек
			for(int i = 0; i < k; i++) {
				Cluster c(my_field.random_point());
				my_search.cl_list.push_back(c);
			}
			// обнуляем кластеры точек
			for (int i = 0; i < my_field.count();i++) {
				my_field.get_point(i).cluster = -1;
			}
			bool changed = true; // поменялась ли принадлежность точки к кластеру
			while (changed == true) {
				changed = fill_clusters(k);
				for (int i = 0; i < k; i++) {
					my_search.cl_list[i].c = my_search.cl_list[i].mCenter();
					my_search.cl_list[i].clear();
				}
			}
			changed = fill_clusters(k);
			my_field.add_search(my_search);
			return my_field.search_count(); // return search_id
		}

		bool print(int search_id, ofstream& print) {
			my_search = my_field.get_search(search_id);
			my_search.print_all(print);
			FILE *f;
			f = popen("gnuplot kmeans.script", "r");
			return pclose(f);
		}
};

class SpanningTreeCl {
	private:
		Field& my_field;
		Cluster_search my_search;
		char** matrix; // граф
		vector <int> a, b; // вспомогательные массивы
	public:
		SpanningTreeCl(Field& f) : my_field(f) {};

		void dfs(int pp, int step, vector <bool>& visit) {
			// уже встречался в пороговой кластеризации
			visit[pp] = true;
			my_search.cl_list[step].add_point(my_field.get_point(pp));
			for (int j = 0; j < my_field.count(); j++) {
				if (visit[j] == false && matrix[pp][j] == 1) {
					dfs(j, step, visit);
				}
			}
		}

		void spanning_tree(ofstream& graph, ofstream& histogram) {
		// Строим дерево, печатаем гистрограмму и дерево
			my_field.fill_matrix(); // заполняем матрицу расстояний
			int size = my_field.count(); // число точек
			// Матрица смежности графа
			matrix = new char*[size]; 
   			for (int i = 0; i < size; i++) {
        		matrix[i] = new char[size];
        		for (int j = 0; j < size;j++)
					matrix[i][j] = 0;
   			}

   			vector <bool> visit(my_field.count(), false); 

			// min_e[i] - минимальное расстояние из вершины i 
			// sel_e[i] - конец наименьшего ребра из i 
			vector<double> min_e (size, 9999);
			vector <int> sel_e (size, -1);
			min_e[0] = 0;
			for (int i = 0; i < size; i++) {
				int v = -1; 
				// Ищем непосещенную вершину, расстояние до которой наименьшее
				// до тех вершин, которые прошли
				for (int j = 0; j < size; ++j)
					if (!visit[j] && (v == -1 || min_e[j] < min_e[v]))
						v = j;
				visit[v] = true;
				if (sel_e[v] != -1) {
					Point p1 = my_field.get_point(v), p2 = my_field.get_point(sel_e[v]);
					histogram << my_field.distances[v][sel_e[v]] << endl;
					graph << p1.getX() << " " << p1.getY()  << endl;
					graph << p2.getX() << " " << p2.getY()  << endl << endl;
					a.push_back(v);
					b.push_back(sel_e[v]);
					matrix[v][sel_e[v]] = matrix[sel_e[v]][v] = 1;
				}
 				// Находим наименьшее ребра для вершины v 
				for (int to = 0; to < size; ++to)
					if (my_field.distances[v][to] < min_e[to]) {
						min_e[to] = my_field.distances[v][to];
						sel_e[to] = v;
					}
			}
			FILE *f;
			f = popen("gnuplot sp_tree.script", "r");
			pclose(f);
		}

		int start(int k) {
			// Нужно построенное дерево
			// k - числов кластеров(количество ребер для удаления)
			my_search.params.push_back("spanning tree");
			my_search.params.push_back(to_string(k));

			vector <bool> visit(my_field.count(), false); 

			// сортируем ребра по длине
			for (int i = 0; i < my_field.count();i++) {
				for (int j = i; j < my_field.count();j++) {
					if (my_field.distances[a[i]][b[i]] < my_field.distances[a[j]][b[j]]) {
						swap(b[i], b[j]);
						swap(a[i], a[j]);
					}
				}
			}

			// Удаляем самые длинные ребра
			for (int i = 0; i < k-1; i++) 
				matrix[a[i]][b[i]] = matrix[b[i]][a[i]] = 0;

			// Выделяем связные компоненты(кластеры) в дереве
			int step = 0, pp = 0;
			while (pp < my_field.count()) { // пока не посетили все точки
				Cluster c(my_field.random_point());
				my_search.cl_list.push_back(c);
				dfs(pp, step, visit);
				for (int i = 0; i <= my_field.count(); i++) {
					if (visit[i] == false) {
						pp = i;
						break;
					}
				}
				step++;
			}

			my_field.add_search(my_search);
			return my_field.search_count(); // return search_id
		}

		bool print(int search_id, ofstream& print) {
			my_search = my_field.get_search(search_id);
			my_search.print_all(print);
			FILE *f;
			f = popen("gnuplot kmeans.script", "r");
			return pclose(f);
		}
};

class Hierarchical_cl {
	private:
		Field& my_field;
		Cluster_search my_search;
		Node *root; // дерево дендрограммы
		int num = 1; // для печати
	public:
		Hierarchical_cl(Field& f) : my_field(f) {};

		void DendrogramBuilder(int p1_id, int p2_id) {
			// обьединяем деревья кластеров p1_id и p2_id
			Node *tree = new Node;
			tree->l = my_field.get_point(p1_id).Tree;
		    tree->r = my_field.get_point(p2_id).Tree;
			tree->aY = my_field.distances[p1_id][p2_id];
			my_field.get_point(p1_id).Tree = tree;
			root = tree;
		}

		int start(double porog) {
			my_search.params.push_back("hierarchical");
			my_search.params.push_back(to_string(porog));
			my_field.fill_matrix();
			int cluster_count = my_field.count();
			int numeration = 0;
			// обьявляем каждую точку кластером
			for (int i = 0; i < my_field.count();i++) {
				my_field.get_point(i).cluster = numeration;
				my_field.get_point(i).Tree = NULL;
				numeration++;
			}

			while (cluster_count > 1) {
				double minimal = 99999;
				int a1 = 0, a2 = 0;
				// Ищем 2 ближайших кластера a1 и a2
				for (int i = 0; i < my_field.count(); i++)  
					for (int l = 0; l < my_field.count(); l++)
							if (my_field.distances[i][l] < minimal && my_field.distances[i][l] > 0) {
								minimal = my_field.distances[i][l];
								a1 = i; a2 = l;
							}
				DendrogramBuilder(a1, a2);
				if (minimal < porog) {
					// обьединяем кластеры a1 и a2
					for (int j = 0; j < my_field.count();j++)
						if (my_field.get_point(j).cluster == a2)
							my_field.get_point(j).cluster = a1;

					for (int i = 0; i < my_field.count(); i++) {
						my_field.distances[a1][i] = my_field.distances[i][a1] = min(my_field.distances[a1][i], my_field.distances[a2][i]);
						my_field.distances[a2][i] = my_field.distances[i][a2] = 99999; // поглощенные кластеры
					}
				}
				cluster_count--;
			}
			// собираем кластеры
			vector <bool> visit(my_field.count(), false);

			for (int i = 0; i < my_field.count(); i++)
				if (visit[i] == false) {
					Point p = my_field.get_point(i);
					Cluster c(p);
					visit[i] = true;
					for (int j = i + 1; j < my_field.count(); j++) {
						if (p.cluster == my_field.get_point(j).cluster) {
							visit[j] = true;
							c.add_point(my_field.get_point(j));
						}
					}
					my_search.cl_list.push_back(c);
				}

			my_field.add_search(my_search);
			return my_field.search_count(); // return search_id	
		}

		void printDendrogram(Node *Tree, ofstream& dendrogram) {
			if (Tree->l != NULL && Tree->r != NULL) {
				printDendrogram(Tree->l, dendrogram);
				printDendrogram(Tree->r, dendrogram);
				dendrogram << (Tree->l)->aX << " " << (Tree->l)->aY << endl;
				dendrogram << (Tree->l)->aX << " " << Tree->aY << endl << endl;
				dendrogram << (Tree->r)->aX << " " << (Tree->r)->aY << endl;
				dendrogram << (Tree->r)->aX << " " << Tree->aY << endl << endl;
				dendrogram << (Tree->l)->aX << " " << Tree->aY << endl;
				dendrogram << (Tree->r)->aX << " " << Tree->aY << endl << endl;
				Tree->aX = (double)((Tree->r)->aX + (Tree->l)->aX)/2;
			}
			else if (Tree->l != NULL){
				printDendrogram(Tree->l, dendrogram);
				dendrogram  << (Tree->l)->aX << " " << (Tree->l)->aY << endl;
				dendrogram << (Tree->l)->aX << " " << Tree->aY << endl << endl;
				dendrogram << num << " " << 0 << endl;
				dendrogram << num << " " << Tree->aY << endl << endl;
				dendrogram << (Tree->l)->aX << " " << Tree->aY << endl;
				dendrogram << num << " " << Tree->aY << endl << endl;
				Tree->aX = (double)((Tree->l)->aX+num)/2;
				num++;
			}
			else if (Tree->r != NULL){
				printDendrogram(Tree->r, dendrogram);
				dendrogram  << num << " " << 0 << endl;
				dendrogram << num << " " << Tree->aY << endl << endl;
				dendrogram << (Tree->r)->aX << " " << (Tree->r)->aY << endl;
				dendrogram << (Tree->r)->aX << " " << Tree->aY << endl << endl;
				dendrogram << (Tree->r)->aX << " " << Tree->aY << endl;
				dendrogram << num << " " << Tree->aY << endl << endl;
				Tree->aX = (double)((Tree->r)->aX+num)/2;
				num++;
			}
			else {
				dendrogram << num << " " << 0 << endl;
				dendrogram << num << " " << Tree->aY << endl << endl;
				dendrogram << num+1 << " " << 0 << endl;
				dendrogram << num+1 << " " << Tree->aY << endl << endl;
				dendrogram << num << " " << Tree->aY << endl;
				dendrogram << num+1 << " " << Tree->aY << endl << endl;
				Tree->aX = (double)(num+num+1)/2;
				num+=2;
			}
		}
		bool print(int search_id, ofstream& print, ofstream& dendrogram) {
			my_search = my_field.get_search(search_id);
			my_search.print_all(print);
			printDendrogram(root, dendrogram); // work only with start(...)
			FILE *f;
			f = popen("gnuplot hierarchical.script", "r");
			return pclose(f);
		}

};
class Forel {
	private:
		Field& my_field;
		Cluster_search my_search;
	public:
		Forel(Field& f) : my_field(f) {};
		int start(double R) {
			// R - радиус кластера
			my_search.params.push_back("FOREL");
			my_search.params.push_back(to_string(R));
			int clusterized = 0; 
			// отмечаем точки, как невыбранные
			for (int i = 0; i < my_field.count();i++) {
				my_field.get_point(i).cluster = -1;
			}

			// Пока есть некластеризованные точки, выделяем кластеры
			while (clusterized < my_field.count()) {
				Point p(0,0); // случайная не кластеризованная точка
				Point center(1,1); // Центр будущего кластера
				Cluster my_cluster(p);
				// Выбираем первую некластеризованную точку p
				for (int i = 0; i < my_field.count();i++) {
					if (my_field.get_point(i).cluster == -1) {
						p = my_field.get_point(i);
						my_cluster.add_point(my_field.get_point(i));
						break;
					}
				}

				// Пока центр кластера меняется, пересчитываем его
				while(fabs(p.getX()-center.getX()) > eps || fabs(p.getY()-center.getY()) > eps) {
					p = my_cluster.mCenter();
					my_cluster.clear();
					for (int i = 0; i < my_field.count();i++) {
						Point new_point = my_field.get_point(i);
						if (p.distanceTo(new_point) < R && new_point.cluster == -1)
							my_cluster.add_point(new_point);
					}

					center = my_cluster.mCenter();
				}
				my_cluster.clear();
				for (int i = 0; i < my_field.count();i++) {
					Point temp_point = my_field.get_point(i);
					if (p.distanceTo(temp_point) < R && temp_point.cluster == -1) {
						my_cluster.add_point(temp_point);
						clusterized++;
						my_field.get_point(i).cluster = 1;

					}
				}
				my_search.cl_list.push_back(my_cluster);
			}
			my_field.add_search(my_search);
			return my_field.search_count(); // return search_id
		}

		bool print(int search_id, ofstream& pointPrint, ofstream& circlePrint) {
			my_search = my_field.get_search(search_id);
			my_search.print_all(pointPrint);
			for (int i = 0; i < my_search.clusterCount(); i++) {
				Point clusterCenter = my_search.get_cluster(i).mCenter();
				circlePrint << clusterCenter.getX() << " " << clusterCenter.getY() << " "<< my_search.getParameter(1) << endl;
			}

			FILE *f;
			f = popen("gnuplot forel.script", "r");
			return pclose(f);
		}
};

class DBSCAN {
	private:
		Field& my_field;
		Cluster_search my_search;
	public:
		DBSCAN(Field& f) : my_field(f) {};

		int start(double R, int minPts) {
			// R - радиус окрестности точки
			// minPts порог для шума
			my_search.params.push_back("DBSCAN");
			my_search.params.push_back(to_string(R));
			my_search.params.push_back(to_string(minPts));

			int cl_num = 0; // номер кластера
			for (int i = 0; i < my_field.count();i++) {
				my_field.get_point(i).cluster = -1;
			}

			for (int i = 0; i < my_field.count();i++) {
				if (my_field.get_point(i).cluster != -1) continue;

				vector <int> Neighbours; // соседи
				for (int j = 0; j < my_field.count();j++) {
					if (my_field.get_point(i).distanceTo(my_field.get_point(j)) < R) // Если лежит в окрестноти 
						Neighbours.push_back(j);
	            }

				if ((int)Neighbours.size() < minPts) { // Проверяем плотность
					my_field.get_point(i).cluster = -2;
					continue;
				}

				cl_num++;
				my_field.get_point(i).cluster = cl_num;
				vector <int> S; // 
				S = Neighbours;

				for (int y = 0; y < (int)S.size(); y++) {
					if (my_field.get_point(S[y]).cluster == -2)
						my_field.get_point(S[y]).cluster = cl_num;
					if (my_field.get_point(S[y]).cluster != -1) continue;

					my_field.get_point(S[y]).cluster = cl_num;
					Neighbours.clear();
					for (int j = 0; j < my_field.count();j++)
						if (my_field.get_point(S[y]).distanceTo(my_field.get_point(j)) < R)
							Neighbours.push_back(j);
					
					if ((int)Neighbours.size() >= minPts)
						for(auto v : Neighbours)
							S.push_back(v);

				}
			}
			// Собираем кластеры
			vector <bool> visit(my_field.count(), false);

			for (int i = 0; i < my_field.count(); i++) {
				if (visit[i] == false) {
					Point p = my_field.get_point(i);
					Cluster c(p);
					visit[i] = true;
					for (int j = i + 1; j < my_field.count(); j++) {
						if (p.cluster == my_field.get_point(j).cluster) {
							visit[j] = true;
							c.add_point(my_field.get_point(j));
						}
					}
					my_search.cl_list.push_back(c);
				}
			}
			my_field.add_search(my_search);
			return my_field.search_count(); // return search_id
		}

		bool print(int search_id, ofstream& print) {
			my_search = my_field.get_search(search_id);
			my_search.print_all(print);
			FILE *f;
			f = popen("gnuplot kmeans.script", "r");
			return pclose(f);
		}
};

class EM_alghoritm {
	private:
		Field& my_field; // поле, которое кластеризуем
		Cluster_search my_search;
		ofstream& animation; // plotting gif
	public:
		EM_alghoritm(Field& f, ofstream& an) : my_field(f), animation(an) {};

		void makeProbabilityClusters(int k, vector<vector<double> >& X, vector <int>& createdClusters) {
			// k - cluster count
			// X - matrix of cluster membership probabilities
			// createdClusters - the order in which clusters go in cl_list

			int it = 0; //iterator
			for (int i = 0; i < my_field.count(); i++) {
				double maxChance = 0; 
				int pointCluster = 0; // cluster with the biggest chance
				for (int j = 0; j < k; j++) {
					if (X[it][j] > maxChance) {
						maxChance = X[it][j];
						pointCluster = j;
					}
				}
				if (createdClusters[pointCluster] > -1) { // if cluster exist
					int cl_id = createdClusters[pointCluster]; // Adding point to cluster
					my_search.cl_list[cl_id].add_point(my_field.get_point(i));
				}
				else { // Creating new cluster
					Cluster c(my_field.get_point(i));
					c.add_point(my_field.get_point(i));
					my_search.cl_list.push_back(c);
					createdClusters[pointCluster] = my_search.cl_list.size()-1;
				}
				it++;
			}
		}
		
		int start(int k, double delta, int Q, bool needToAnimate) {
			my_search.params.push_back("EM");
			my_search.params.push_back(to_string(k));
			my_search.params.push_back(to_string(delta));

			vector <int> createdClusters(k, -1); // the order in which clusters go in cl_list

			vector <vector <double> > C(2); // matrix (q*k) of expected values 
			vector <double> R(2*k, 1); // k diagonal covariance matrixes (2*2)
			vector <double> W(k, 0); // weights (k*1)
			// size - point count
			int size = my_field.count();
			vector<vector<double> > X(size, vector<double>(k, 0)); // output - our probabilities (size * k)

			// Giving random probabilities
			for (int i = 0; i < 2; i++) {
				C[i].resize(k);
				for (int j = 0; j < k; j++) {
					C[i][j] = (double)(rand() % 1000)/1000;
				}
			}

			for (int i = 0; i < k; i++) {
				W[i] = (double)1/k;
			}

			double llh = 0, new_llh = 2*delta;
			// Beginning
			int step = 0;
			while (step < Q && abs(new_llh-llh) > delta) {
				step++;
				// E step
				llh = new_llh;
				new_llh = 0;
				vector <vector <double> > tC(2, vector<double>(k, 0)); // temp matrixes
				vector <double> tR(2*k, 0);
				vector <double> tW(k, 0);
				for (int i = 0; i < size; i++) {
					double sumpi = 0;
					Point pi = my_field.get_point(i);
					for (int j = 0; j < k; j++) {
						double sigma = ((pi.getX()-C[0][j])*(pi.getX()-C[0][j])/R[2*j]) + (pi.getY()-C[1][j]) * (pi.getY()-C[1][j])/R[2*j+1];
						X[i][j] = W[j] * exp(-sigma/2)/sqrt(2*PI)/sqrt(R[2*j] * R[2*j+1]);
						sumpi+=X[i][j];
					}
					new_llh+= log(sumpi);
					for (int j = 0; j < k; j++) {
						X[i][j]/=sumpi;
						tW[j]+=X[i][j];
						tC[0][j]+= pi.getX() * X[i][j]; 
						tC[1][j]+= pi.getY() * X[i][j]; 
					}
				}
				// M step
				for (int j = 0; j < k; j++) {
					for (int l = 0; l < 2; l++) {
						C[l][j] = tC[l][j]/tW[j];
					}
					for (int i = 0; i < size; i++) {
						Point pi = my_field.get_point(i);
						tR[2*j] += (pi.getX()-C[0][j])*(pi.getX()-C[0][j])*X[i][j];
						tR[2*j+1] += (pi.getY()-C[1][j]) * (pi.getY()-C[1][j])*X[i][j];
					}
					for (int i = 0; i < k; i++) {
						W[i] = tW[i]/size;
					}
					R[2*j] = tR[2*j]/tW[j];
					R[2*j+1] = tR[2*j+1]/tW[j];
				}

				///// ANIMATION
				if (needToAnimate) {
					if (step % 10 == 0) { // often animation causes segmentation fault
						makeProbabilityClusters(k, X, createdClusters);
						animate(k);
						my_search.deleteClusters();
					}
				}
				/////
			}
			// End, we have got X matrix
			cout << step << " steps, deviation = " << abs(new_llh - llh) << endl;
			makeProbabilityClusters(k, X, createdClusters); // based on matrix X, we make clusters
			
			if (needToAnimate) {
				animate(k);
				FILE *f;
				f = popen("gnuplot animation.script", "r");
				pclose(f);
			}

			my_field.add_search(my_search);
			return my_field.search_count(); // return search_id
		}
		void animate(int k) {
			for (int i = 0; i < k; i++) {
				Matrix2D covClusterMatrix;
				vector <double> eighen;
				Point center = my_search.cl_list[i].mCenter();
				covClusterMatrix.make_covariance(my_search.cl_list[i].points, center);
				eighen = covClusterMatrix.eighensystem();
				animation << endl << center.getX() << " " << center.getY() << " "
				<< 2*sqrt(eighen[4]) << " " << 2*sqrt(eighen[5]) << " " << acos(eighen[0]/sqrt(eighen[0]*eighen[0]+eighen[1]*eighen[1]))*180/PI << endl;
			}
			animation << endl << endl;
			my_search.print_all(animation);
			animation << endl;
		}

		bool print(int search_id, ofstream& print, ofstream& ellipses) {
			// Печатаем кластеры, матрицу ковариации для каждого кластера
			// и эллипс
			my_search = my_field.get_search(search_id);
			my_search.print_all(print);

			int k = stoi(my_search.getParameter(1)); // число кластеров
			for (int i = 0; i < k; i++) {
				Matrix2D covClusterMatrix;
				vector <double> eighen;
				Point center = my_search.cl_list[i].mCenter();
				covClusterMatrix.make_covariance(my_search.cl_list[i].points, center);
				eighen = covClusterMatrix.eighensystem();
				ellipses << center.getX() << " " << center.getY() << " " << covClusterMatrix.get(0,0) << " " <<
				covClusterMatrix.get(0, 1) << " " << covClusterMatrix.get(1,0) << " " << covClusterMatrix.get(1,1) << " ";
				ellipses << 2*sqrt(eighen[4]) << " " << 2*sqrt(eighen[5]) << " " << acos(eighen[0]/sqrt(eighen[0]*eighen[0]+eighen[1]*eighen[1]))*180/PI << endl;
			}

			FILE *f;
			f = popen("gnuplot em.script", "r");
			return pclose(f);
		}
};

class Control {
	// управляющий класс
	private:
		Field f;
		int rootId = -1;
		stringstream stringOut;
	public:
		void print_group(int i) { 
			// вывести все точки группы
			ofstream gnu("data/group.txt");
			f.gr_list[i].print_group(gnu, false);
			stringOut << "recorded to 'data/group.txt'" << endl;
		}
		void print_all() {
			f.print_all();
			stringOut << "recorded to 'data/field.txt'" << endl;
		}
		void clear() {
			f.clear();
		}
		void archieve() {
			f.archieve();
		}
		void group_list() {
			f.group_list();
		}
		void search_list() {
			f.search_list();
		}
		void genrnd(int min, int max) {
			// равномерная генерация N точек с координатами от min до max
			if (f.getMode()) {
				stringOut << "Error. Clustering mode" << endl;
				return;
			}
			Group gr(0, 0, 0, 0, 0);
			for (int i = 0; i < 1000; i++) {
				Point p(rand() % (max-min) + min, rand() % (max-min) + min);
				gr.add_point(p);
			}
			f.add(gr);
			stringOut << "group: generated | id = " << f.gr_list.size() - 1 << endl;
		}
		void move(int id, double x, double y) {
			if (f.getMode()) { // Проверка состояния поля
				stringOut << "Error. Clustering mode" << endl;
				return;
			}
			f.gr_list[id].move(x, y);
		}
		void rotate(int id, double alfa) {
			if (f.getMode()) { // Проверка состояния поля
				stringOut << "Error. Clustering mode" << endl;
				return;
			}
			f.gr_list[id].rotate(alfa);
		}
		void scale(int id, double r1, double r2) {
			if (f.getMode()) { // Проверка состояния поля
				stringOut << "Error. Clustering mode" << endl;
				return;
			}
			f.gr_list[id].scale(r1, r2);
		}
		void create_group(double x, double y, double disp1, double disp2, int N, int client) {
			// Генерируем группу и добавляем ее в поле

			if (f.getMode()) { // Проверка состояния поля
				stringOut << "Error. Clustering mode" << endl;
				return;
			}

			Group gr(x, y, disp1, disp2, N);
			f.add(gr);

			stringOut << "group: created | id = " << f.gr_list.size() - 1 << endl;
		}
		void kmeans(int k, int search_id) {
			// либо кластеризует и печатает (search_id == -1)
			// либо только печатает (search_id > 0)
			if (!f.getMode()) {
				stringOut << "Error. Generation mode" << endl;
				return;
			}
			Kmeans clust(f);
			if (search_id == -1) {
				int id = clust.start(k);
				if (id >= 0) {
					stringOut << "search: added | id = " << id << endl;
				}
			}
			ofstream print("data/clusters.txt");
			if (!clust.print(search_id, print)) {
				stringOut << "image created" << endl;
			}

		}
		void porog(double d, int search_id) {
			if (!f.getMode()) {
				stringOut << "Error. Generation mode" << endl;
				return;
			}
			Porog clust(f);
			if (search_id == -1) {
				int id = clust.start(d);
				if (id >= 0) {
					stringOut << "search: added | id = " << id << endl;
				}
			}
			ofstream print("data/clusters.txt");
			if (!clust.print(search_id, print)) {
				stringOut << "image created" << endl;
			}
		}
		void spanning_tree(int search_id) {
			// ! Здесь ведется диалог с клиентом
			// ! Без кластеризации дерево не печатает
			int k; // число кластеров
			if (!f.getMode()) {
				stringOut << "Error. Generation mode" << endl;
				return;
			}
			SpanningTreeCl clust(f);
			if (search_id == -1) {
				ofstream graph("data/graph.txt");
				ofstream histogram("data/histogram.txt");
				clust.spanning_tree(graph, histogram);
				stringOut << "Введите число кластеров:" << endl;
				cin >> k;
				int id = clust.start(k);
				if (id >= 0) {
					stringOut << "search: added | id = " << id << endl;
				}
			}
			ofstream print("data/clusters.txt");
			if (!clust.print(search_id, print)) {
				stringOut << "image created" << endl;
			}
		}
		void hierarchical_cl(double porog, int search_id) {
			if (!f.getMode()) {
				stringOut << "Error. Generation mode" << endl;
				return;
			}
			Hierarchical_cl clust(f);
			if (search_id == -1) {
				int id = clust.start(porog);
				if (id >= 0) {
					stringOut << "search: added | id = " << id << endl;
				}
			}
			ofstream print("data/clusters.txt");
			ofstream dendrogram("data/dendrogram.txt");
			if (!clust.print(search_id, print, dendrogram)) {
				stringOut << "image created" << endl;
			}
		}
		void forel(double R, int search_id) {
			if (!f.getMode()) {
				stringOut << "Error. Generation mode" << endl;
				return;
			}
			Forel clust(f);
			if (search_id == -1) { // нужно кластеризировать
				int id = clust.start(R);
				if (id >= 0) {
					stringOut << "search: added | id = " << id << endl;
				}
			}
			ofstream print("data/clusters.txt");
			ofstream circles("data/circles.txt");
			if (!clust.print(search_id, print, circles)) {
				stringOut << "image created" << endl;
			}

		}
		void dbscan (double R, int minPts, int search_id) {
			if (!f.getMode()) {
				stringOut << "Error. Generation mode" << endl;
				return;
			}
			DBSCAN clust(f);
			if (search_id == -1) {
				int id = clust.start(R, minPts);
				if (id >= 0) {
					stringOut << "search: added | id = " << id << endl;
				}
			}
			ofstream print("data/clusters.txt");
			if (!clust.print(search_id, print)) {
				stringOut << "image created" << endl;
			}
		}
		void emAlghoritm(int k, double delta, int Q, int search_id, bool needToAnimate) {
			if (!f.getMode()) {
				stringOut << "Error. Generation mode" << endl;
				return;
			}
			ofstream animation("data/animation.txt");
			EM_alghoritm clust(f, animation);
			if (search_id == -1) {
				int id = clust.start(k, delta, Q, needToAnimate);
				if (id >= 0) {
					stringOut << "search: added | id = " << id << endl;
				}
			}
			ofstream print("data/clusters.txt");
			ofstream ellipses("data/ellipses.txt");
			if (!clust.print(search_id, print, ellipses)) {
				stringOut << "image created" << endl;
			}
		}
		void print_clusters(int id) {
			//f.searches[id].print_all();
			stringOut << "recorded to 'data/clusters.txt'" << endl;
		}
		void print_pca(int id) {
			//f.searches[id].print_pca();
			stringOut << "recorded to 'data/pca.txt'" << endl;
		}
		void load(string filename) {
			ifstream file(filename);
			f.load(file);
		}
		void changeMode(int client) {
			if (client == rootId) {
				f.changeMode();
				stringOut << "changed" << endl;
			}
			else {
				stringOut << "No permission" << endl;
			}
		}
        void identify(int client, string command) {
        	LOG << "FROM [" << client << "]: " << command << endl;
        	if (rootId < 0)
        		rootId = client;
        	vector <string> tokens;
        	int words = 0;
            string b;
			size_t pos = 0;
			while ((pos = command.find(" ")) != string::npos) {
   				 tokens.push_back(command.substr(0, pos));
   				 words++;
   				 command.erase(0, pos + 1);
			}
			tokens.push_back(command);
			words++;
            if (tokens[0] == "genrnd" && words > 2)
                genrnd(stoi(tokens[1]), stoi(tokens[2]));
            else if (tokens[0] == "group" && words > 5) 
            	create_group(stod(tokens[1]), stod(tokens[2]), stod(tokens[3]), stod(tokens[4]), stoi(tokens[5]), client);
            else if (tokens[0] == "print"&& words > 1)
            	print_group(stoi(tokens[1]));
            else if (tokens[0] == "gr")
            	print_all();
             else if (tokens[0] == "mode")
            	changeMode(client);	
            else if (tokens[0] == "archieve")
            	archieve();
            else if (tokens[0] == "clear")
            	clear();
            else if (tokens[0] == "cl" && words > 1)
            	print_clusters(stoi(tokens[1]));
            else if (tokens[0] == "pca" && words > 1)
            	print_pca(stoi(tokens[1]));
            else if (tokens[0] == "gr_list")
            	group_list();
            else if (tokens[0] == "searches")
            	search_list();
            else if (tokens[0] == "move" && words > 3)
            	move(stoi(tokens[1]), stod(tokens[2]), stod(tokens[3]));
            else if (tokens[0] == "rotate")
            	rotate(stoi(tokens[1]), stod(tokens[2]));
            else if (tokens[0] == "kmeans" && words > 1)
            	kmeans(stoi(tokens[1]), -1);
            else if (tokens[0] == "porog" && words > 1)
            	porog(stod(tokens[1]), -1);
            else if (tokens[0] == "sp_tree")
            	spanning_tree(-1);
            else if (tokens[0] == "hcl" && words > 1) {
				hierarchical_cl(stod(tokens[1]), -1);
            }
            else if (tokens[0] == "forel" && words > 1)
				forel(stod(tokens[1]), -1);
			else if (tokens[0] == "dbscan" && words > 1)
				dbscan(stod(tokens[1]), stoi(tokens[2]), -1);
			else if (tokens[0] == "em") {
				if (words > 3)
            		emAlghoritm(stoi(tokens[1]), stod(tokens[2]), stoi(tokens[3]), -1, true);
            	else 
            		emAlghoritm(stoi(tokens[1]), 0.0001, 300, -1, false);
			}
            else if (tokens[0] == "load") {
            	if (words == 1) load("data/ARCHIEVE.txt");
            	else load(tokens[1]);
            }
            else if (tokens[0] == "help")  {
                ifstream help ("help.txt");
                while (!help.eof()) {
                    getline(help, b);
                    stringOut << b << endl;
                }
            }
            else if (tokens[0] != "exit") {
                stringOut << tokens[0] << " command not found or incorrect parameters." << endl << "Use 'help' to see commands." << endl;
            }
            writeToClient(client, stringOut);
        }
        
};

class Server {
	private:
		string ipAddress;
		int port;
		Control c;

		int createSocket() {
			int sock = socket (AF_INET, SOCK_STREAM, 0);
			if (sock < 0) {
				cerr << "Server: cannot create socket" << endl;
        		return -1;
			}
			struct sockaddr_in addr;
			// Заполняем адресную структуру и
		    // связываем сокет с любым адресом
		    addr.sin_family = AF_INET;
		    addr.sin_port = htons(port);
		    addr.sin_addr.s_addr = htonl(INADDR_ANY);

		    int err = bind(sock,(struct sockaddr*)&addr,sizeof(addr));
		    if (err < 0) {
		    	cerr << "Server: cannot bind socket" << endl;
        		return -1;
		    }
		    // Создаем очередь на 3 входящих запроса соединения
		    err = listen(sock,3);
		    if (err < 0) {
		        cerr << "Server: listen queue failure" << endl;
		        return -1;
		    }
		    return sock;
		}

	public:
		Server(string ip, int prt) : ipAddress(ip), port(prt) {}

		void run() {
			char buf[4096];
			sockaddr_in  client;
			socklen_t  size;

			int sock = createSocket();
			if (sock < 0) {
		    	return;
		    }
			// Подготавливаем множества дескрипторов каналов ввода-вывода.
		    // Для простоты не вычисляем максимальное значение дескриптора,
		    // а далее будем проверять все дескрипторы вплоть до максимально
		    // возможного значения FD_SETSIZE.
		    pollfd  act_set[100];
		    act_set[0].fd = sock;
		    act_set[0].events = POLLIN;
		    act_set[0].revents = 0;
		    int num_set = 1;

	    	while(true) {
	    		// Проверим, не появились ли данные в каком-либо сокете.
		        // В нашем варианте ждем до фактического появления данных.
		        int ret = poll (act_set, num_set, -1);
		        if (ret < 0) {
		            cerr << "Server: poll  failure" << endl;
		            break;
		        }
		        if (ret > 0) {
		        	for (int i = 0; i < num_set; i++) {
		        		if (act_set[i].revents & POLLIN) {
		        			cout << "get POLLIN at fd" << act_set[i].fd << endl;
		        			act_set[i].revents &= ~POLLIN;
		        			// Новое подключение
		        			if (i == 0) {
		        				int new_sock = accept(act_set[i].fd,(struct sockaddr*)&client,&size);
		        				cout << "new client at port" <<  ntohs(client.sin_port) << endl;
		        				if (num_set < 100) {
		        					act_set[num_set].fd = new_sock;
		        					act_set[num_set].events = POLLIN;
		        					act_set[num_set].revents = 0;
		        					num_set++;
		        				} 
		        				else {
		        					cout << "no more sockets for client" << endl;
		        					close(new_sock);
		        				}
		        			}
		        			else {
		                   		// пришли данные в уже существующем соединени
		                   		int nbytes = recv(act_set[i].fd, buf, 65536, 0);
		                   		c.identify(act_set[i].fd, string(buf));
								   cout << string(buf) << endl;
		                        if (nbytes <= 0) { // ошибка или конец данных
		                      		cout << "get stop" << endl;
		                      		close (act_set[i].fd);
		                      		if (i < num_set-1) {
		                      			act_set[i] = act_set[num_set - 1];
		                      			num_set--;
		                      			i--;
		                      		}
		                      	}
							}  
						}
	           		}
    			}
	    	}
		}
};	

int main() {
	Server my_server("127.0.0.1", 5555);
	my_server.run();
    return 0;
}
