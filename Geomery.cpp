#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <map>
#include <iomanip>
#include <fstream>
#include <queue>
#include <stack>
#include <iomanip>
#include <numeric>
#include <limits>

using namespace std;

typedef long long ll;

#define PI 3.14159265358979323846 

const int MOD = (int)1e9 + 7;

class Point {
public:

	Point(double long x = 0, double long y = 0, double long z = 0, bool isPoint = true) {
		this->x = x, this->y = y, this->z = z;
		if (isPoint)this->toPoint();
		else		this->toVector();
	}

	/*Instance of objects in this class*/
	bool isInstancePoint() {
		return "Point" == this->instance;
	}
	bool isInstanceVector() {
		return "Vector" == this->instance;
	}
	void toPoint() {
		this->instance = "Point";
	}
	void toVector() {
		this->instance = "Vector";
	}

	/*Addresses*/
	double long operator[](int iterator){
		if (iterator == 0) {
			return this->x;
		}
		else if (iterator == 1) {
			return this->y;
		}
		return this->z;
	}
	double long operator[](string item) {
		if (item == "x" || item == "X") {
			return this->x;
		}
		else if (item == "y" || item == "Y") {
			return this->y;
		}
		return this->z;
	}
	double long operator[](char item) {
		if (item == 'x' || item == 'X') {
			return this->x;
		}
		else if (item == 'y' || item == 'Y') {
			return this->y;
		}
		return this->z;
	}

	/*
	Operations
	Operator *  -dot
	Operator %  -cross3d
	Operator &  -cross2d
	Operator +  -add
	Operator -  -minus
	Operator ^  -multiply (double long)
	Operator =  -assign
	Operator == -eq
	*/
	//Operators
	double long operator*(Point item) { // -dot
		return this->x * item['X'] + this->y * item['Y'] + this->z * item['Z'];
	}
	double long operator&(Point item) { // cross 2d
		return item['X'] * this->y - this->x * item['Y'];
	}

	Point operator%(Point item) { //-cross 3d
		return { this->y * item['Z'] - this->z * item['Y'], this->x * item['Z'] - this->z * item['X'], this->x * item['Y'] - this->y * item['X'] };
	}
	void operator%=(Point item) { //-cross 3d
		double long tempX = this->y * item['Z'] - this->z * item['Y'], tempY = this->x * item['Z'] - this->z * item['X'], tempZ = this->x * item['Y'] - this->y * item['X'];
		this->x = tempX, this->y = tempY, this->z = tempZ;
	}

	Point operator+(Point item) { //-sum
		return { this->x + item['X'], this->y + item['Y'], this->z + item['Z'] };
	}
	void operator+=(Point item) { //-sum
		this->x += item['X'], this->y += item['Y'], this->z += item['Z'];
	}

	Point operator-(Point item) { //-dif
		return { this->x - item['X'], this->y - item['Y'], this->z - item['Z'] };
	}
	void operator-=(Point item) { //-dif
		this->x -= item['X'], this->y -= item['Y'], this->z -= item['Z'];
	}
	
	Point operator^(double long amount) { //-mlpy const
		return { this->x * amount, this->y * amount, this->z * amount };
	}
	void operator^=(double long amount) { //-mlpy const
		this->x *= amount, this->y *= amount, this->z *= amount;
	}

	bool operator==(Point item) { //-isEqual
		return (this->x == item['X'] && this->y == item['Y'] && this->z == item['Z']);
	}
	void operator=(Point item) { //-assign
		this->x = item['X'], this->y = item['Y'], this->z = item['Z'];
	}

	//stdio stream
	friend ostream& operator<<(ostream& output, const Point& pt) {
		output << pt.x << " " << pt.y;
		return output;
	}
	friend istream& operator>>(istream& input, Point& pt) {
		input >> pt.x >> pt.y;
		return input;
	}

	//Methods
	void assign(double long x = 0, double long y = 0, double long z = 0, bool isPoint = true) {
		this->x = x, this->y = y, this->z = z; 
		if (isPoint)this->toPoint();
		else		this->toVector();
	}
	void rotateRight(double long angle) {
		this->rotateLeft(-angle);
	}
	void rotateLeft(double long angle) {
		double long cs = cos(angle * PI / 180), sn = sin(angle * PI / 180);
	

		double long tempX = this->x * cs - this->y * sn;
		double long tempY = this->x * sn + this->y * cs;
		this->x = tempX, this->y = tempY;
	}

	double long norm() {
		return sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
	}

	double long sin2D(Point item) {
		return this->pCross2D(item) / (this->norm() * item.norm());
	}
	double long cos2D(Point item) {
		return (this->x * item['X'] + this->y * item['Y'] + this->z * item['Z']) / (this->norm() * item.norm());
	}

	vector<double long>getCords3D(){
		return {this->x, this->y, this->z};
	}
	vector<double long>getCords2D() {
		return { this->x, this->y };
	}
	pair<double long, double long>getByPairCords2D() {
		return { this->x, this->y };
	}


private:
	double long x = 0, y = 0, z = 0;
	string instance = "Point";

	double long pCross2D(Point item) { return item['X'] * this->y - this->x * item['Y']; }

};

class Line {
public:

	Line(Point A, Point B) { // Points init
		Point AB = B - A;
		this->A = AB['Y'];
		this->B = -AB['X'];
		this->C = -(this->A * A['X'] + this->B * A['Y']);
		ll GCD = gcd(gcd(this->A, this->B), this->C);
		this->A /= GCD, this->B /= GCD, this->C /= GCD;
		
	}
	Line(double long A = 0, double long B = 0, double long C = 0) {
		this->A = A, this->B = B, this->C = C;
		ll GCD = gcd(gcd(this->A, this->B), this->C);
		this->A /= GCD, this->B /= GCD, this->C /= GCD;
	}

	Point intersection(Line item) {
		double long P = item['A'], Q = item['B'], T = item['C'], A = this->A, B = this->B, C = this->C;
		return {((C*Q-B*T)/(B*P-A*Q)),((C*P-A*T)/(A*Q-B*P))};
	}

	Line orthogonal(Point item) {
		return {-this->B, this->A, -(-this->B * item['X'] + this->A * item['Y']) };
	}

	Line bisector(Line with, Point LeftBorder, Point RightBorder, bool externalAngle = false) {
		double long P = this->A, Q = this->B, T = this->C, R = with['A'], U = with['B'], K = with['C'];
		Point N1(P, Q), N2(R, U);
		double long n1 = N1.norm(), n2 = N2.norm();

		double long a = n2*P - n1*R, b = n2*Q - U*n1, c = T*n2 - K*n1;
		double long a_ = n2 * P + n1 * R, b_ = n2 * Q + U * n1, c_ = T * n2 + K * n1;
		Line pBisector(a, b, c), pBisector2(a_,b_,c_), Border(LeftBorder, RightBorder);

		Point borderPoint = pBisector.intersection(Border);

		Point LR = LeftBorder - RightBorder;
		Point LbP = LeftBorder - borderPoint;
		Point RbP = RightBorder - borderPoint;

		if (abs(LR.norm() - (LbP.norm() + RbP.norm())) <= 0.0000001 ) {
			if (!externalAngle) return pBisector;
			else				return pBisector2;
		}
		else {
			if (!externalAngle) return pBisector2;
			else				return pBisector;
		}
	}

	Point GetNormal() {
		return { this->A, this->B };
	}

	void inverse() {
		this->A = -this->A;
		this->B = -this->B;
		this->C = -this->C;
	}
	
	bool isCollinearTo(Line item) {
		return (item['A'] == this->A && item['B'] == this->B || -item['A'] == this->A && -item['B'] == this->B);
	}
	bool isOnLine(Point item) {
		return (abs(item['X'] * this->A + item['Y'] * this->B + this->C) <= 0.000001);
	}

	/*Addresses*/
	double long operator[](int iterator) {
		if (iterator == 0) {
			return this->A;
		}
		else if (iterator == 1) {
			return this->B;
		}
		return this->C;
	}
	double long operator[](string item) {
		if (item == "a" || item == "A") {
			return this->A;
		}
		else if (item == "b" || item == "B") {
			return this->C;
		}
		return this->C;
	}
	double long operator[](char item) {
		if (item == 'a' || item == 'A') {
			return this->A;
		}
		else if (item == 'b' || item == 'B') {
			return this->B;
		}
		return this->C;
	}

	/*Operators*/
	void operator=(Line item) {
		this->A = item['A'], this->B = item['B'], this->C = item['C'];
	}
	bool operator==(Line item) {
		return (this->A == item['A'] && this->B == item['B'] && this->C == item['C'] || this->A == -item['A'] && this->B == -item['B'] && this->C == -item['C']);
	}
	//stdio stream

	friend ostream& operator<<(ostream& output, const Line& ln) {
		output << '(' << ln.A << ")x + (" << ln.B << ")y + (" << ln.C << ") = 0";
		return output;
	}
	friend istream& operator>>(istream& input, Line& ln) {
		input >> ln.A >> ln.B >> ln.C;
		return input;
	}

private:
	double long A = 0, B = 0, C = 0;

	ll gcd(ll a, ll b)
	{
		return b ? gcd(b, a % b) : a;
	}
};

class Cycle {
public:

	Cycle(double long x = 0, double long y = 0, double long r = 0) {
		this->x = x, this->y = y, this->r = r;
	}
	Cycle(Point A, Point B, Point C, bool isInTriangle = false) {
		if (isInTriangle) {
			Line AB(A, B), BC(B, C), AC(A, C);

			Line fBisector = AB.bisector(BC, A, C);
			Line sBisector = BC.bisector(AC, A, B);

			Point center = fBisector.intersection(sBisector);

			Point N = AB.intersection(AB.orthogonal(center));

			this->x = center['X'], this->y = center['Y'], this->r = (N - center).norm();
			return;
		}
		Line AB(A, B), BC(B, C);
		Point Mab = (A + B) ^ 0.5;
		Point Mbc = (B + C) ^ 0.5;

		Line TMab = AB.orthogonal(Mab);
		Line TMbc = AB.orthogonal(Mbc);

		Point Center = TMab.intersection(TMbc);

		this->x = Center['X'];
		this->y = Center['Y'];
		this->r = (B - Center).norm();
		
	}
	Cycle(Point A, double long r = 0) {
		this->x = A['X'], this->y = A['Y'], this->r = r;
	}

	/* Methods */
	double long chord(Point A) {
		double long dist = 0;
		Point C(this->x, this->y);
		dist = (C - A).norm();
		return sqrt((this->r * this->r - dist*dist) * 4);
	}

	vector<Point> LineIntersection(Line item) {
		double long A = item['A'], B = item['B'], C = item['C'], R = this->r;

		Point Center(this->x, this->y);
		Line Ort = item.orthogonal(Center);
		Point CI = item.intersection(Ort);
		double long MR = (CI - Center).norm();
		if (R >= MR && R - MR < 0.000001) {
			return {CI};
		}
		else if (R < MR) {
			return {};
		}
		Point Vnorm = Ort.GetNormal();
		double long t = this->chord(CI)/(2*Vnorm.norm());

		return {(Vnorm ^ (t)) + CI, (Vnorm ^ (-t)) + CI};
	}

	vector<Point> CycleIntersection(Cycle item) {
		
		if (this->p_eq(item)) {
			return { {0,0}, {0,0}, {0,0} };// Inf Points;
		}

		double long P = this->x, Q = this->y, R = this->r;
		double long K = item['X'], T = item['Y'], U = item['Z'];
		double long A = 2 * (P-K);
		double long B = -2 * (T-Q);
		double long C = K*K - P*P - Q*Q + R*R + T*T - U*U;
		Line INTCC(A, B, C);

		return this->LineIntersection(INTCC);
	}

	bool onCycle(Point item) {
		double long calc = ((this->x - item['X']) * (this->x - item['X']) + (this->y - item['Y']) * (this->y - item['Y'])) - (this->r * this->r);
		return calc > 0 && calc < 0.000001;
	}

	/*Addresses*/
	double long operator[](int iterator) {
		if (iterator == 0) {
			return this->x;
		}
		else if (iterator == 1) {
			return this->y;
		}
		return this->r;
	}
	double long operator[](string item) {
		if (item == "x" || item == "X") {
			return this->x;
		}
		else if (item == "y" || item == "Y") {
			return this->y;
		}
		return this->r;
	}
	double long operator[](char item) {
		if (item == 'x' || item == 'X') {
			return this->x;
		}
		else if (item == 'y' || item == 'Y') {
			return this->y;
		}
		return this->r;
	}

	//stdio stream
	friend istream& operator>>(istream& input, Cycle& cle) {
		input >> cle.x >> cle.y >> cle.r;
		return input;
	}

	//operators
	bool operator==(Cycle item) {
		return (item['X'] == this->x && item['Y'] == this->y && item['R'] == this->r);
	}


private:
	double long x = 0, y = 0, r = 0;

	bool p_eq(Cycle item) {
		return (item['X'] == this->x && item['Y'] == this->y && item['R'] == this->r);
	}
};
