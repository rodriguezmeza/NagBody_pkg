/*
 * buildmap.c: generate and compile client transformation program.
 */

#include "stdinc.h"
#include "datatypes.h"
#include "getparam.h"
#include "buildmap.h"
#include <ctype.h>
#include <unistd.h>

local void codemap(stream sstr, string *names, string *exprs, string *types,
		   string tmap);

local string mapdefs[][2];
local string expdefs[][2];

//  buildmap: generate extendbody, computemap, and computetime routines;
//  compile and link with snapmap main program.
//  ____________________________________________________________________

void buildmap(string prog, string *names, string *exprs, string *types,
	      string tmap, string prec, int ndim, bool cleanup)
{
  string src, cmd;
  stream sstr;

  asprintf(&src, "%s.c", prog);
  sstr = stropen(src, "w");
  codemap(sstr, names, exprs, types, tmap);
  fclose(sstr);
  asprintf(&cmd, "%s %s %s -D%s -DNDIM=%d -o %s %s.c %s/lib/snapmap_%c%d.o "
	   "-lNBody -lClib -lgsl -lgslcblas -lm", getenv("ZCC"), 
	   getenv("ZCCFLAGS"), getenv("ZLDFLAGS"), prec, ndim, prog, prog,
	   getenv("ZENOPATH"), tolower(prec[0]), ndim);
  eprintf("[%s: compiling snapmap: %s]\n", getprog(), cmd);
  if (system(cmd) != 0)
    error("%s: command \"%s\" failed\n", getprog(), cmd);
  if (cleanup && unlink(src) != 0)
    error("%s: can't unlink %s\n", getprog(), src);
  free(src);
  free(cmd);
}

//  getmapdefs: return a pointer to the mapdefs table.
//  __________________________________________________

string *getmapdefs(void)
{
  return ((string *) mapdefs);
}

//  codemap: generate actual code for mapping operation.
//  ____________________________________________________

local void codemap(stream sstr, string *names, string *exprs, string *types,
		   string tmap)
{
  string zss = getenv("ZENO_SAFE_SELECT"), *tp, *np, *ep, tn, *arrexp;
  int i;

  fprintf(sstr, "#include \"stdinc.h\"\n");
  fprintf(sstr, "#include \"mathfns.h\"\n");
  fprintf(sstr, "#include \"vectdefs.h\"\n");
  if (zss == NULL || strne(zss, "FALSE"))	// check unless set to FALSE
    fprintf(sstr, "#define SafeSelect TRUE\n");
  fprintf(sstr, "#include \"phatbody.h\"\n\n");
  for (i = 0; mapdefs[i][0] != NULL; i++)
    fprintf(sstr, "#define %-8s %s(_p)\n", mapdefs[i][0], mapdefs[i][1]);
  fprintf(sstr, "\n");
  for (i = 0; expdefs[i][0] != NULL; i++)
    fprintf(sstr, "#define %-8s %s\n", expdefs[i][0], expdefs[i][1]);
  fprintf(sstr, "\n");
  if (types != NULL) {
    for (tp = types, np = names; *tp != NULL; tp++, np++) {
      if (*np == NULL)
	error("%s.buildmap: more types than names\n", getprog());
      tn = type_name(*tp);
      fprintf(sstr, "#define %s(b)  Select%c%s"
	      "(b, phatbody[NewBodyFields+%d].offset)\n",
	      *np, toupper(tn[0]), tn + 1, (int) (tp - types));
    }
    fprintf(sstr, "\n");
    fprintf(sstr, "void extendbody(void)\n");
    fprintf(sstr, "{\n");
    for (tp = types, np = names; *tp != NULL; tp++, np++)
      fprintf(sstr, "  new_field(&phatbody[NewBodyFields+%d], "
	      "\"%s\", \"%s\");\n", (int) (tp - types), *tp, *np);
    fprintf(sstr, "  new_field(&phatbody[NewBodyFields+%d],"
	    " NULL, NULL);\n", (int) (tp - types));
    fprintf(sstr, "}\n\n");
  } else
    fprintf(sstr, "void extendbody(void)\n{ }\n\n");
  fprintf(sstr,
	  "void computemap(bodyptr _q, bodyptr _p, real t, int i, int n)\n");
  fprintf(sstr, "{\n");
  for (ep = exprs, np = names; *ep != NULL; ep++, np++) {
    if (*np == NULL)
      error("%s.buildmap: more exprs than names\n", getprog());
    if (index(*ep, ';') == NULL)
      fprintf(sstr, "  %s(_q) = (%s);\n", *np, *ep);
    else {
      arrexp = burststring(*ep, ";");
      for (i = 0; arrexp[i] != NULL; i++)
	fprintf(sstr, "%s(_q)[%d] = (%s);\n", *np, i, arrexp[i]);
    }
  }
  fprintf(sstr, "}\n\n");
  fprintf(sstr, "real computetime(real t, int n)\n");
  fprintf(sstr, "{\n");
  fprintf(sstr, "  return (%s);\n", tmap != NULL ? tmap : "t");
  fprintf(sstr, "}\n");
}

//  mapdefs: mapping between identifiers used in expressions and macro names.
//  _________________________________________________________________________

local string mapdefs[][2] = {
  { "pos",     "Pos"     },
  { "x",       "PosX"    },
  { "y",       "PosY"    },
  { "z",       "PosZ"    },
  { "vel",     "Vel"     },
  { "vx",      "VelX"    },
  { "vy",      "VelY"    },
  { "vz",      "VelZ"    },
  { "m",       "Mass"    },
  { "phi",     "Phi"     },
  { "acc",     "Acc"     },
  { "ax",      "AccX"    },
  { "ay",      "AccY"    },
  { "az",      "AccZ"    },
  { "smooth",  "Smooth"  },
  { "rho",     "Rho"     },
  { "entf",    "EntFunc" },
  { "uint",    "Uintern" },
  { "udot",    "UdotInt" },
  { "udotrad", "UdotRad" },
  { "udotvis", "UdotVis" },
  { "tau",     "Tau"     },
  { "type",    "Type"    },
  { "birth",   "Birth"   },
  { "death",   "Death"   },
  { "key",     "Key"     },
  { "keyarr",  "KeyArr"  },
  { "aux",     "Aux"     },
  { "auxv",    "AuxVec"  },
  { "auxvx",   "AuxVecX" },
  { "auxvy",   "AuxVecY" },
  { "auxvz",   "AuxVecZ" },
  { "auxarr",  "AuxArr"  },
  { NULL,      NULL      }
};

//  expdefs: quantities derived from basic variables.
//  _________________________________________________

local string expdefs[][2] = {
  { "r",       "absv(pos)"                        },
  { "R",       "rsqrt(x*x + y*y)"                 },
  { "v",       "absv(vel)"                        },
  { "vr",      "(dotvp(pos,vel)/r)"               },
  { "vt",      "rsqrt(dotvp(vel,vel) - rsqr(vr))" },
  { "etot",    "(phi + 0.5*dotvp(vel,vel))"       },
  { "jx",      "(vy*z - vz*y)"                    },
  { "jy",      "(vz*x - vx*z)"                    },
  { "jz",      "(vx*y - vy*x)"                    },
  { "jtot",    "rsqrt(jx*jx + jy*jy + jz*jz)"     },
  { NULL,       NULL                              }
};

#ifdef TESTBED

string defv[] = {		";Invoke buildmap routine",
  "prog=sm",			";Name of program to compile",
  "names=XPlot,YPlot,Color",	";Names of output variables",
  "exprs=x;y;1+i%2",		";Expressions for output variables.",
				";Bound variables, depending on input, are:",
				  SNAPMAP_BODY_VARS ".",
  "types=" FloatType "," FloatType "," IntType,
				";Data types of output vars",
  "tmap=t+1",			";Expression used to map time",
				";Bound variables are: "
				  SNAPMAP_TIME_VARS ".",
  "prec=SINGLEPREC",		";Specify precision option",
  "ndim=3",			";Number of space dimensions",
  "VERSION=1.5",		";Josh Barnes  2 February 2015",
  NULL,
};

int main(int argc, string argv[])
{
  string *names, *exprs, *types;

  initparam(argv, defv);
  names = burststring(getparam("names"), ",");
  exprs = burststring(getparam("exprs"), ";");
  if (! strnull(getparam("types")))
    types = burststring(getparam("types"), ",");
  else
    types = NULL;
  buildmap(getparam("prog"), names, exprs, types,
	   getparam("tmap"), getparam("prec"), getiparam("ndim"), FALSE);
  return (0);
}

#endif
