#include "stdlib.h"
#include "math.h"
#include <stdio.h>
#include <iostream>
#include <limits>
#include <sqlite3.h> 
#include <sstream>
#include <string>

#define BRDF_SAMPLING_RES_THETA_H       90
#define BRDF_SAMPLING_RES_THETA_D       90
#define BRDF_SAMPLING_RES_PHI_D         360

const char *brdffile = "/Users/niverik2k/Downloads/polyurethane-foam.binary";
const char *outputfile = "/Users/niverik2k/Downloads/polyurethane-foam.binary";
double* brdf;

static int callback(void *NotUsed, int argc, char **argv, char **azColName)
{
   int i;
   for(i=0; i <argc; i++)
   {
      printf("%s = %s\n", azColName[i],argv[i] ? argv[i] : "NULL");
   }
printf("\n");
return 0;
}


// I use this function when I want to create a table, I will have to incorporate it into the main function eventually.
int createtable()
{
   sqlite3 *db;
   char *zErrMsg = 0;
   int rc;
   char *sql;

   rc = sqlite3_open("test.db", &db);

   if( rc )
   {
      fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
      return(0);
   }
   else
   {
      fprintf(stderr, "Opened database successfully\n");
   }

   /* Create SQL statement */
   sql= const_cast<char *>("CREATE TABLE BRDFDATA(" \
      "ID INT PRIMARY	KEY	NOT NULL," \
      "VALUE		TEXT	NOT NULL);");
   /*Execute SQL statement */
   rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);
   if(rc != SQLITE_OK)
   {
      fprintf(stderr, "SQL error: %s\n", zErrMsg);
      sqlite3_free(zErrMsg);
   }
   else
   {
      fprintf(stdout, "Table created sucessfully\n");
   }
   sqlite3_close(db);
   return 0;
}

// Insert data item into BRDFDATA table, it takes individual records with an ID and VALUE.
void RunInsertParamSQL(sqlite3 *db, int brdfid, double brdfvalue)
{
  if (!db)
    return;

  char *zErrMsg = 0;
  sqlite3_stmt *stmt;
  const char *pzTest;
  char *szSQL;

  //This is where the sqlite statement is generated, the ? marks allow me to bind variables later.
  szSQL = const_cast<char *>("insert into BRDFDATA (ID, VALUE) values (?,?)");
  
  //I prepare the statement, to get it in the correct format.
  int rc = sqlite3_prepare(db, szSQL, strlen(szSQL), &stmt, &pzTest);
  
  // Check to make sure everything went through ok with the prep, befor binding and commiting.
  if( rc == SQLITE_OK ) {
    // bind the value 
    sqlite3_bind_int(stmt, 1, brdfid);
    sqlite3_bind_double(stmt, 2, brdfvalue);

    // commit 
    sqlite3_step(stmt);
    sqlite3_finalize(stmt);
  }
}
   

// Read BRDF data, I open the brdf file, and assign it to the brdf double array pointer that is n units large.
// I return n to generate an index for the values later
int read_brdf(const char *brdffile, double* &brdf)
{
	FILE *f = fopen(brdffile, "rb");
	if (!f)
		return -2;

	int dims[3];
	fread(dims, sizeof(int), 3, f);
	int n = dims[0] * dims[1] * dims[2];
	if (n != BRDF_SAMPLING_RES_THETA_H *
		 BRDF_SAMPLING_RES_THETA_D *
		 BRDF_SAMPLING_RES_PHI_D / 2) 
	{
		fprintf(stderr, "Dimensions don't match\n");
		fclose(f);
		return -2;
	}

	brdf = (double*) malloc (sizeof(double)*3*n);
	fread(brdf, sizeof(double), 3*n, f);
	fclose(f);
	return n;
}

// I take the brdf pointer and sqlite3 pointer to write the brdf data to the sqlite database.
int writefile(double* &brdf, sqlite3 *db)//, const char *outputfile)
{
   int n = read_brdf(brdffile, brdf);
   int ctr;

   for(ctr=0;ctr <=n;ctr++)
   {
      RunInsertParamSQL(db, ctr, brdf[ctr]);
   }
   printf("finished");
   return 0;
}


//int n = 20;
//int main()
//{
//	writefile(brdf);
//
//	return 0;
//}

// Opens the sqlite databse, then uses the writefile() function to write to the sqlite database.
int main()
{
   sqlite3 *db;
   char *zErrMsg = 0;
   int rc;
   char *sql;

   rc = sqlite3_open("test.db", &db);

   if( rc )
   {
      fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
      return(0);
   }
   else
   {
      fprintf(stderr, "Opened database successfully\n");
   }
   //The commented code takes a long time to generate, only uncomment it if you intend to write to the database.
   //writefile(brdf, db);

   sqlite3_close(db);
   return 0;
}

