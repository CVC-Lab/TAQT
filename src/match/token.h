#ifndef _DX2RAW_TOK_H
#define _DX2RAW_TOK_H

#ifdef __cplusplus
extern "C" {
#endif

#define ID 1  
#define NUM 2 
#define COMMENT 3

#define MAXLEN 256

extern char yylval[MAXLEN];

int parse_comment(); 

int parse_id(); 

int parse_fnum(); 

#ifdef __cplusplus
}
#endif

#endif

