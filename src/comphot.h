#ifndef COMPHOT_H
#define COMPHOT_H
#pragma once

typedef struct {
	const char* offsetimage;
	const char* fixedimage;
	const char* flatimage;
	int cenx;
	int ceny;
	float apradius;
} ComphotConfig;

void process( const ComphotConfig* config );

#endif // COMPHOT_H

