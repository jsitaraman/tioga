#ifndef LINKLIST_H
#define LINKLIST_H
//
// This file is part of the Tioga software library
//
// Tioga  is a tool for overset grid assembly on parallel distributed systems
// Copyright (C) 2015 Jay Sitaraman
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

#include "codetypes.h"

void deallocateLinkList(DONORLIST *temp);
void deallocateLinkList2(INTEGERLIST *temp);
void deallocateLinkList3(INTEGERLIST2 *temp);
void deallocateLinkList4(INTERPLIST2 *temp);
void insertInList(DONORLIST **donorList,DONORLIST *temp1);

#endif // LINKLIST_H
