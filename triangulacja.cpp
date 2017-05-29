#include <iostream>
#include <cmath>
#include <time.h>
#include <float.h>
#include <fstream>
#include <GL/glut.h>
#include <stdlib.h>
#include <limits>
using namespace std; 

struct TRIANGLE
{
    int p1, p2, p3;
};

struct EDGE
{
    int p1, p2;
};

struct XYZ
{
    double x, y, z;
};

////////////////////////////////////////////////////////////////////////
// CircumCircle() :
//   Zwraca true jeśli punkt o współrzędnych (xp,yp) znajduje się w środku (x1,y1), (x2,y2), (x3,y3)
//   Środek okręgu opisanego jest zwracany w postaci zmiennych(xc,yc) and the radius r
//   Punkt leżący na okręgu liczy się jako leżący w okręgu.
////////////////////////////////////////////////////////////////////////

int CircumCircle(double xp, double yp, double x1, double y1, double x2, double y2, double x3, double y3, double &xc, double &yc, double &r)
{
    double m1, m2, mx1, mx2, my1, my2;
    double dx, dy, rsqr, drsqr;
    double mindif = DBL_EPSILON;

//wynajdywanie punktów o jednakowym położeniu
    if(abs(y1 - y2) < mindif && abs(y2 - y3) < mindif)
    return(false);
    if(abs(y2-y1) < mindif)
    { 
        m2 = - (x3 - x2) / (y3 - y2);
        mx2 = (x2 + x3) / 2.0;
        my2 = (y2 + y3) / 2.0;
        xc = (x2 + x1) / 2.0;
        yc = m2 * (xc - mx2) + my2;
    }
    else if(abs(y3 - y2) < mindif)
    { 
        m1 = - (x2 - x1) / (y2 - y1);
        mx1 = (x1 + x2) / 2.0;
        my1 = (y1 + y2) / 2.0;
        xc = (x3 + x2) / 2.0;
        yc = m1 * (xc - mx1) + my1;
    }
    else
    {
         m1 = - (x2 - x1) / (y2 - y1); 
         m2 = - (x3 - x2) / (y3 - y2); 
         mx1 = (x1 + x2) / 2.0; 
         mx2 = (x2 + x3) / 2.0;
         my1 = (y1 + y2) / 2.0;
         my2 = (y2 + y3) / 2.0;
         xc = (m1 * mx1 - m2 * mx2 + my2 - my1) / (m1 - m2); 
         yc = m1 * (xc - mx1) + my1; 
    }
    dx = x2 - xc;
    dy = y2 - yc;
    rsqr = dx * dx + dy * dy;
    r = sqrt(rsqr); 
    dx = xp - xc;
    dy = yp - yc;
    drsqr = dx * dx + dy * dy;
    return((drsqr <= rsqr) ? true : false);
}
///////////////////////////////////////////////////////////////////////////////
// Triangulate() :
//   Przebieg triangulacji.
//   nv punktów zostaje wprowadzone do kontenera 'pxyz'.
//   Zostaje zwrócona lista 'ntri' trójkątnych powierzchni w kontenerze 'v'.
//   Trójkąty muszą zostać uporządkowane zgodnie z ruchem wskazówek zegara.
//   Kontenerowi zawierającemu trójkąty 'v' powinia zostać przydzielona pamięć na pomieszczenie '3 * nv' obiekrów.
//   Kontener zawierający punkty 'pxyz' musi pomieścić jeszcze punkty tworzące trójkąt zawierający wszystkie punkty.
//   Kontener zawierający punkty musi zostać posortowany wraz ze wzrastającą wartością x. z wykorzystaniem funkcji qsort(p,nv,sizeof(XYZ),XYZCompare);
///////////////////////////////////////////////////////////////////////////////

int Triangulate(int nv, XYZ pxyz[], TRIANGLE v[], int &ntri)
{
    int *complete = NULL;
    EDGE *edges = NULL; 
    EDGE *p_EdgeTemp;
    int nedge = 0;
    int trimax, emax = 200;
    int status = 0;
    int inside;
    int i, j, k;
    double xp, yp, x1, y1, x2, y2, x3, y3, xc, yc, r;
    double xmin, xmax, ymin, ymax, xmid, ymid;
    double dx, dy, dmax; 

//Przypisanie pamięci dla powstających trójkątów przy pomocy operatora new, rezerwującego ilość pamięci równą (wielkość zmiennej) * (ilość zmiennych)
    trimax = 4 * nv;
    complete = new int[trimax];
//Analogiczne przypisanie pamięci dla krawędzi
    edges = new EDGE[emax];

// Znajdź punkty skrajne.
// Potrzebne do tworzenia trójkąta zawierającego wszystkie punkty.

    xmin = pxyz[0].x;
    ymin = pxyz[0].y;
    xmax = xmin;
    ymax = ymin;
    for(i = 1; i < nv; i++)
    {
        if (pxyz[i].x < xmin) xmin = pxyz[i].x;
        if (pxyz[i].x > xmax) xmax = pxyz[i].x;
        if (pxyz[i].y < ymin) ymin = pxyz[i].y;
        if (pxyz[i].y > ymax) ymax = pxyz[i].y;
    }
    dx = xmax - xmin;
    dy = ymax - ymin;
    dmax = (dx > dy) ? dx : dy;
    xmid = (xmax + xmin) / 2.0;
    ymid = (ymax + ymin) / 2.0;

//Tworzenie trójkąta zawierającego wszystkie wprowadzone punkty.
//Znajduje się on na pierwszej pozycji w liście utworzonych trójkątów.

    pxyz[nv+0].x = xmid - 20 * dmax;
    pxyz[nv+0].y = ymid - dmax;
    pxyz[nv+1].x = xmid;
    pxyz[nv+1].y = ymid + 20 * dmax;
    pxyz[nv+2].x = xmid + 20 * dmax;
    pxyz[nv+2].y = ymid - dmax;
    v[0].p1 = nv;
    v[0].p2 = nv+1;
    v[0].p3 = nv+2;
    complete[0] = false;
    ntri = 1;

//Wprowadzanie kolejnych punktów jeden po drugim do 

    for(i = 0; i < nv; i++)
    {
        xp = pxyz[i].x;
        yp = pxyz[i].y;
        nedge = 0;

//Tworzenie krawędzi:
// Jeśli punkt (xp,yp) znajduje się wewnątrz trójkąta, to jego trzy krawędzie zostają dodane do zbioru, a trójkąt zostaje usunięty.


        for(j = 0; j < ntri; j++)
        {
            if(complete[j])
            continue;
            x1 = pxyz[v[j].p1].x;
            y1 = pxyz[v[j].p1].y;
            x2 = pxyz[v[j].p2].x;
            y2 = pxyz[v[j].p2].y;
            x3 = pxyz[v[j].p3].x;
            y3 = pxyz[v[j].p3].y;
            inside = CircumCircle(xp, yp, x1, y1, x2, y2, x3, y3, xc, yc, r);
            if (xc + r + DBL_MIN < xp)
                complete[j] = true;
            if(inside)
            {
//Sprawdznie, czy nie przekroczono maksymalnego rozmiaru pamięci przydzielonej na zbiór krawędzi, jeśli tak to powiększenie jej rozmiaru.
                if(nedge + 3 >= emax)
                {
                    emax += 100;
                    p_EdgeTemp = new EDGE[emax];
                    for (int i = 0; i < nedge; i++) 
                    { 
                        p_EdgeTemp[i] = edges[i];   
                    }
                delete []edges;
                edges = p_EdgeTemp;
                }
                edges[nedge+0].p1 = v[j].p1;
                edges[nedge+0].p2 = v[j].p2;
                edges[nedge+1].p1 = v[j].p2;
                edges[nedge+1].p2 = v[j].p3;
                edges[nedge+2].p1 = v[j].p3;
                edges[nedge+2].p2 = v[j].p1;
                nedge += 3;
                v[j] = v[ntri-1];
                complete[j] = complete[ntri-1];
                ntri--;
                j--;
            }
        }

        for(j = 0; j < nedge - 1; j++)
        {
            for(k = j + 1; k < nedge; k++)
            {
                if((edges[j].p1 == edges[k].p2) && (edges[j].p2 == edges[k].p1))
                {
                    edges[j].p1 = -1;
                    edges[j].p2 = -1;
                    edges[k].p1 = -1;
                    edges[k].p2 = -1;
                }

            }
        }

        for(j = 0; j < nedge; j++) 
        {
            if(edges[j].p1 < 0 || edges[j].p2 < 0)
            continue;
            v[ntri].p1 = edges[j].p1;
            v[ntri].p2 = edges[j].p2;
            v[ntri].p3 = i;
            complete[ntri] = false;
            ntri++;
        }
    }

//Usuwanie punktów z trójkąta zawierającego wszystkie punkty.
//Są to punkty dla których wartość i jest większa niż ilość wierzchołków zdefiniowanych nv.

    for(i = 0; i < ntri; i++)
    {
        if(v[i].p1 >= nv || v[i].p2 >= nv || v[i].p3 >= nv)
        {
            v[i] = v[ntri-1];
            ntri--;
            i--;
        }
    }
    delete[] edges;
    delete[] complete;
    return 0;
} 

int XYZCompare(const void *v1, const void *v2)
{
    XYZ *p1, *p2;
    
    p1 = (XYZ*)v1;
    p2 = (XYZ*)v2;
    if(p1->x < p2->x)
        return(-1);
    else if(p1->x > p2->x)
        return(1);
    else
        return(0);
}

int random(int n)//konieczne bo rand() jest typu int a trzeba wprowadzić do double
{
    return rand()%n; 
}

int main( int argc, char * argv[] )
{
    bool check;
    int n_MaxPoints;
    cout << "Podaj maksymalną ilość punktów." << endl;
    do
    {
        cin >> n_MaxPoints;
        check = cin.good();
        if(check == false)
        {
            cin.clear();
            cin.ignore( 1000, '\n' );
            cout << "Podaj liczbę." << endl;
        }
    }while (check == false);
    TRIANGLE *v = NULL;
    int max = 10;//maksymalny rozmiar pamięci rezerwowanej za jednym razem
    XYZ *p = new XYZ[max];//rezerwowanie pamięci operatorem new 
    XYZ *p_Temp = NULL;   
    int nv = 0;
    int X, Y;
    int i;
    int ntri = 0;
    double x, y, z;
    bool b_Ok = false;
    srand(time(NULL));
    nv = 0;
    p = new XYZ[max];
    while (nv != n_MaxPoints)
    {
        do//generowanie punktów
        {
            b_Ok = true;
            x = (double)random(1001);
            y = (double)random(1001);
            for(int n_Cpt = 0; n_Cpt <= nv; n_Cpt++)
            {
                if((x == p[n_Cpt].x) && (y == p[n_Cpt].y)) b_Ok = false;
            }//zabezpieczenie przed utworzeniem jednakowych punktów
        }while(!b_Ok);
        if (nv >= max)
        {
            max = max * 2;//podwojenie rozmiaru 
            p_Temp = new XYZ[max]; 
            for (int i = 0; i < nv; i++)
            {
                p_Temp[i] = p[i];  
            }
            delete []p;  
            p = p_Temp; 
        }   
        p[nv].x = x * 1.0;
        p[nv].y = y * 1.0;
        nv++;
    }
    p_Temp = new XYZ[nv + 3]; 
    for (int i = 0; i < nv; i++)
    {
        p_Temp[i] = p[i];      
    }
    delete []p;           
    p = p_Temp;
    v = new TRIANGLE[3 * nv]; 
    qsort(p, nv, sizeof(XYZ), XYZCompare);
    Triangulate(nv, p, v, ntri);
    fstream vertices;
    vertices.open("vertices.txt");
    glutInit( & argc, argv );
    glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGB );
    glutInitWindowSize( 1000, 1000 );
    glutCreateWindow( "Triangulacja" );
    glClearColor( 0.0, 0.0, 0.0, 0.0 );
    glClear( GL_COLOR_BUFFER_BIT );
    glColor3d( 0.0, 1.0, 0.0 );
    float xp, yp, xk, yk;
    for(int i = 0; i < ntri; i++)
    {
        cout << ((float)p[v[i].p1].x-500)/500 << " " << ((float)p[v[i].p1].y-500)/500 << " "<< ((float)p[v[i].p2].x-500)/500 << " " << ((float)p[v[i].p2].y-500)/500 << endl;
        vertices << ((float)p[v[i].p1].x-500)/500 << " " << ((float)p[v[i].p1].y-500)/500 << " " << ((float)p[v[i].p2].x-500)/500 << " " << ((float)p[v[i].p2].y-500)/500 << endl;
        xp = ((float)p[v[i].p1].x-500)/500;
        yp = ((float)p[v[i].p1].y-500)/500;
        xk = ((float)p[v[i].p2].x-500)/500;
        yk = ((float)p[v[i].p2].y-500)/500;
        glBegin( GL_LINES );
        glVertex2d( xp , yp );
        glVertex2d( xk , yk );
        glEnd();
        cout << ((float)p[v[i].p2].x-500)/500 << " "<< ((float)p[v[i].p2].y-500)/500 << " " << ((float)p[v[i].p3].x-500)/500 << " " << ((float)p[v[i].p3].y-500)/500 << endl;
        vertices << ((float)p[v[i].p2].x-500)/500 << " " << ((float)p[v[i].p2].y-500)/500 << " " << ((float)p[v[i].p3].x-500)/500 << " " << ((float)p[v[i].p3].y-500)/500 << endl;
        xp = ((float)p[v[i].p2].x-500)/500;
        yp = ((float)p[v[i].p2].y-500)/500;
        xk = ((float)p[v[i].p3].x-500)/500;
        yk = ((float)p[v[i].p3].y-500)/500;
        glBegin( GL_LINES );
        glVertex2d( xp , yp );
        glVertex2d( xk , yk );
        glEnd();
        cout << ((float)p[v[i].p3].x-500)/500 << " " << ((float)p[v[i].p3].y-500)/500 << " " << ((float)p[v[i].p1].x-500)/500 << " " << ((float)p[v[i].p1].y-500)/500 << endl;
        vertices << ((float)p[v[i].p3].x-500)/500 << " " << ((float)p[v[i].p3].y-500)/500 << " " << ((float)p[v[i].p1].x-500)/500 << " " << ((float)p[v[i].p1].y-500)/500 << endl;
        xp = ((float)p[v[i].p3].x-500)/500;
        yp = ((float)p[v[i].p3].y-500)/500;
        xk = ((float)p[v[i].p1].x-500)/500;
        yk = ((float)p[v[i].p1].y-500)/500;
        glBegin( GL_LINES );
        glVertex2d( xp , yp );
        glVertex2d( xk , yk );
        glEnd();
    }
    glFlush();
    glutSwapBuffers();
    glutMainLoop();
    vertices.close();
    delete []p;
    delete []v;
    p = NULL;
    v = NULL;
    return 0;
}
