/*  scamp_protocol.h */

/*
 * Copyright (c) 2024 Daniel Marks

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
 */

#ifndef __SCAMP_PROTOCOL_H
#define __SCAMP_PROTOCOL_H

#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <string>

class SCAMP_protocol {
private:
	unsigned char val;
	unsigned char table[256];
	char ss[3];
public:
	SCAMP_protocol() { 
	}
	~SCAMP_protocol() {};
	void init();
};

#endif  /* _SCAMP_PROTOCOL_H */