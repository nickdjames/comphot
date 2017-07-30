#ifndef COMPHOT_H
#define COMPHOT_H
#pragma once

// Windows specific export definition
#ifdef _WIN32
#define LIB_API __declspec(dllexport)
#else
#define LIB_API
#endif

typedef LIB_API struct ComphotConfig {
	const char* offsetimage;
	const char* fixedimage;
	const char* flatimage;
	const char *icqtemplate;
	int cenx;
	int ceny;
	int border;
	int photcheck;
	float apradius;
} ComphotConfig;

LIB_API void process( const ComphotConfig* config );

#endif // COMPHOT_H

