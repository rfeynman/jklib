#ifndef __VALU_CLASS_INCL__
#define  __VALU_CLASS_INCL__

const int FORML = 128 + 2;
//const int FakeEvent = 253;

class Value_container
{
public:
  int * formula;
  Value * form_val;

  void fta_zero() const
  {
    for(int i=0; i<256+2; i++) fta[i]=0;
  }
  
  void table_zero() const
  {
    for(int i=0; i<128+2; i++)
      {
	rowval[i]=0;
	rowtim[i]=0;
      }
  }
  
  void formula_zero() const
  {
    for(int i=0; i<FORML; i++) formula[i]=0;
  }
  
  int * cmd_a;
  Value * cmd;
  
  double * rowval;
  Value * val_row;
  
  float  * rowtim;
  Value * val_tim;
  
  int  * table;
  Value * val_table;
  
  int  * fta;
  Value * val_fta;
  
  int   * connEv;
  Value * val_connect;

  Value * val;
  Value * val_a;
  Value * val_b;
  Value * val_c;
  Value * val_f;
  Value * val_g;
  Value * val_h;
  //Value * val_x;
  Value * val_doit;
  Value * val_discA;
  Value * val_Trigg;
  Value * val_int[257];

  Value * val_dzero;
  Value * val_d500th;
  
  Value * val_Yes;
  
  void formAdd(int& i, char c1, char c2, char c3, char c4) const
  {
    formula[2+i++] = c1;
    formula[2+i++] = c2;
    formula[2+i++] = c3;
    formula[2+i++] = c4;
  }
  
  void formConst(int& i, int c) const
  {
    int kk = abs(c);
    char * kkk = (char *) &kk;
    char c1 = kkk[3];
    char c2 = kkk[2];
    char c3 = kkk[1];
    char c4 = c < 0 ? ';' : ':';
    formAdd(i, c1, c2, c3, c4);
  }
  
public:
  Value_container()
  {
    formula = new int[ FORML ];
    //form_int = (int *) formula;
    form_val = new Value( formula, FORML, PRESERVE_ARRAY );
    formula_zero();
    
    cmd_a = new int[4];
    cmd  =  new Value( cmd_a, 4, TAKE_OVER_ARRAY );
    
    rowval = new double[ 128 ];
    val_row = new Value( rowval, 128, TAKE_OVER_ARRAY );
    for(int i = 0; i < 128; i++) rowval[i]=0.;
    
    rowtim = new float[ 128 ];
    val_tim=new Value(rowtim, 128, TAKE_OVER_ARRAY);
    for(int i = 0; i < 128; i++) rowtim[i] = 0.;
    
    table = new int[ 258 ];
    val_table = new Value( table, 258, TAKE_OVER_ARRAY );
    for(int i = 0; i < 258; i++) table[i] = 0;
    
    fta = new int[ 258 ];
    val_fta = new Value( fta, 258, TAKE_OVER_ARRAY );
    for(int i = 0; i < 258; i++) fta[i] = 0;
    
    connEv = new int[ 2 ];
    for(int i =0; i< 2; i++) connEv[i]=0;
    val_connect = new Value( connEv, 2, TAKE_OVER_ARRAY );
    
    
    val_a   = new Value("A");
    val_b   = new Value("B");
    val_c   = new Value("C");
    val_f   = new Value("F");
    val_g   = new Value("G");
    val_h   = new Value("H");
    //val_x   = new Value("X");
    val_doit = new Value("DoIt");
    val_discA = new Value("DisconnectAll");
    val_Trigg = new Value("Trigger");
    val_Yes = new Value("Yes");
    for(int i=0; i<257; i++) val_int[i]   = new Value(i);
    
    val_dzero = new Value(double(0.0));
    val_d500th = new Value(double(0.002));
  }
  ~Value_container()
  {
    Unref(form_val);
    Unref(val_row);
    Unref(val_tim);
    Unref(val_fta);
    Unref(val_table);
    Unref(val_f);
    Unref(val_a);
    Unref(val_b);
    Unref(val_g);
    Unref(val_h);
    //Unref(val_x);
    Unref(val_doit);
    Unref(val_discA);
    Unref(val_Trigg);
    Unref(val_Yes);
    for(int i=0; i<257; i++) Unref(val_int[i]);
    
    Unref(val_dzero);
    Unref(val_d500th);
    Unref(val_connect);
  }
};


extern const Value_container valu;
#endif
