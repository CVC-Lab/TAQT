/* regular definitions */
delim  [ \t\n\,]
ws     {delim}+
letter   [A-Za-z]
digit    [0-9]
id       {letter}({letter}|{digit})*
fnumber  (\-|\.)?{digit}+(\.{digit}+)?((E|e)[+\-]?{digit}+)?(\.)?

%{

#include "token.h"

char yylval[MAXLEN];   
%}

%%
{ws}      { /* no action and no return */ }
{id}      {return(parse_id());}
{fnumber} {return(parse_fnum());}  
"#"       {return(parse_comment());}  
%%

int parse_comment() {
   char c;
   while(1) {
     c = yyinput();
     if(c != '\n');
     else {
       return COMMENT;
     }
   }
}

int parse_id() {
   int i, n;
   n = yyleng;
   if(n >= MAXLEN-1) n = MAXLEN-1;
   for(i = 0; i < n; i++) {
     yylval[i] = yytext[i];
   }
   yylval[i] = '\0';
 
   return ID;
}

int parse_fnum() {
   int i, n;
   n = yyleng;
   if(n >= MAXLEN-1) n = MAXLEN-1;
   for(i = 0; i < n; i++) {
     yylval[i] = yytext[i];
   }
   yylval[i] = '\0';
 
   return NUM;
}
