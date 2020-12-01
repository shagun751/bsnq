#include <iostream>
using namespace std;

class sol
{
public:
	int i;

	sol(int in1){ i = in1; }

	void pr (int j)
	{
		cout << "yes " << i << " " << j << "\n";		
	}
};

typedef void * opqObj;

extern "C"
{
	opqObj GetObject(int);
	void display(opqObj, int);
}

opqObj GetObject(int in1)
{
	sol *s = new sol(in1);
	return (opqObj)s;
}

void display(opqObj foo, int n)
{
	sol *s = (sol *)foo;

	s->pr( n );
}