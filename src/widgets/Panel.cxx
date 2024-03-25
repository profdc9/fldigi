// ----------------------------------------------------------------------------
// Panel.cxx
//
// Copyright (C) 2007-2011
//		Dave Freese, W1HKJ
//
// This file is part of fldigi.
//
// Fldigi is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 or later of the License.
//
// Fldigi is distributed in the hope that it will be useful,
// but WITHOUT ARY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with fldigi.	If not, see <http://www.gnu.org/licenses/>.
// ----------------------------------------------------------------------------

//#include <config.h>

#include <stdlib.h>
#include <stdio.h>

#include <FL/Fl.H>
#include <FL/Enumerations.H>
#include <FL/Fl_Window.H>

#include "Panel.h"
#include <config.h>

// Drag the edges that were initially at oldx,oldy to newx,newy:
// pass -1 as oldx or oldy to disable drag in that direction:
//#define DEVEL_DEBUG 1

int Panel::orgx()
{
	int oldx = w();
	int *p = sizes() + 8;
	for (int i=children(); i--; p += 4)
		if (p[1] < oldx) oldx = p[1];
	return oldx;
}

int Panel::orgy()
{
	int oldy = h();
	int *p = sizes() + 8;
	for (int i=children(); i--; p += 4)
		if (p[3] < oldy) oldy = p[3];
	return oldy;
}

void Panel::position(int oix, int oiy, int newx, int newy) {
#ifdef DEVEL_DEBUG
printf("oix %3d, oiy %3d, nux %3d, nuy %3d\n", oix, oiy, newx, newy);
#endif
	Fl_Widget* const* a = array();
	int *p = sizes();
#ifdef DEVEL_DEBUG
printf("p0 %3d, p1 %3d, p2 %3d, p3 %3d\n", p[0], p[1], p[2], p[3]);
printf("p4 %3d, p5 %3d, p6 %3d, p7 %3d\n", p[0], p[1], p[2], p[3]);
#endif
	p += 8; // skip group & resizable's saved size
	for (int i=children(); i--; p += 4) {
		Fl_Widget* o = *a++;
		if (o == resizable()) continue;
		int X = o->x();
		int Y = o->y();
		int W = o->w();
		int H = o->h();
		int R = X + W;
		int B = Y + H;

		if (o == resizable())
			continue;
		if (oix > -1) {
			int t = p[0];
			if ((t == oix) || (t>oix && X<newx) || (t<oix && X>newx)) X = newx;
			t = p[1];
			if ((t == oix) || (t>oix && R<newx) || (t<oix && R>newx)) R = newx;
		}
		if (oiy > -1) {
			int t = p[2];
			if ((t == oiy) || (t>oiy && Y<newy) || (t<oiy && Y>newy)) Y = newy;
			t = p[3];
			if ((t == oiy) || (t>oiy && B<newy) || (t<oiy && B>newy)) B = newy;
		}
		o->damage_resize(X, Y, R-X, B-Y);
	}
}

// move the lower-right corner (sort of):
void Panel::resize(int X,int Y,int W,int H) {
	// remember how much to move the child widgets:
	int *p = sizes();
	int OX = x();
	int OY = y();
	int OW = w();
	int OH = h();
	int OR = OX + OW;
	int OB = OY + OH;
	float dw = 1.0 * W / OW;
	float dh = 1.0 * H / OH;
	// resize this (skip the Fl_Group resize):
	Fl_Widget::resize(X,Y,W,H);
	// find x, y, w, h of resizable:
	int RX = X + (p[4] - p[0]);
	int RY = Y + (p[6] - p[2]);
	int RR = X + W - (p[1] - p[5]);
	int RB = Y + H - (p[3] - p[7]);
	int NW = RR - RX;
	int NH = RB - RY;
	int R;
	int B;
	int xx;
	int yy;
	// move everything to be on correct side of new resizable:
	Fl_Widget * const *a = array();
	Fl_Widget * o = 0;
	p += 8;
	for (int i=children(); i--;) {
		o = *a++;
		if (o == resizable()) {
			o->resize(RX, RY, NW, NH);
		} else {
			xx = X;
			if (o->x() != OX)
				xx = (o->x() - OX) * dw + X + 0.5;
			if (xx > RR) xx = RR;
			if (o->x() + o->w() == OR) {
				R = X + W;
				if (xx != X && xx < RX) xx = RX;
			} else {
				R = xx + o->w() * dw + 0.5;
				if (R < xx) R = xx;
				if (xx <= RX && R < RX) R = RX;
				if (xx <= RX && R > RR ) R = RR;
			}

			yy = Y;
			if (o->y() != OY)
				yy = (o->y() - OY) * dh + Y + 0.5;
			if (yy > RB) yy = RB;
			if (o->y() + o->h() == OB) {
				B = Y + H;
				if (yy != Y && yy < RY) yy = RY;
			} else {
				B = yy + o->h() * dh + 0.5;
				if (B < yy) B = yy;
				if (yy <= RY && B < RY) B = RY;
				if (yy <= RY && B > RB) B = RB;
			}

			o->resize(xx,yy,R-xx,B-yy);
		}
		p += 4; // next child sizes array
// do *not* call o->redraw() here! If you do, and the tile is inside a 
// scroll, it'll set the damage areas wrong for all children!
	}
}

static void set_cursor(Panel *t, Fl_Cursor c) {
	static Fl_Cursor cursor;
	if (cursor == c || !t->window()) return;
	cursor = c;
	t->window()->cursor(c);
}

static Fl_Cursor cursors[4] = {
	FL_CURSOR_DEFAULT,
	FL_CURSOR_WE,
	FL_CURSOR_NS,
	FL_CURSOR_MOVE};

int Panel::handle(int event) {
	static int sdrag;
	static int sdx, sdy;
	static int sx, sy;
#define DRAGH 1
#define DRAGV 2
#define GRABAREA 4

	int mx = Fl::event_x();
	int my = Fl::event_y();

	switch (event) {

	case FL_MOVE:
	case FL_ENTER:
	case FL_PUSH: {
		int mindx = 100;
		int mindy = 100;
		int oldx = 0;
		int oldy = 0;
		Fl_Widget*const* a = array();
		int *q = sizes();
		int *p = q + 8;
		for (int i=children(); i--; p += 4) {
			Fl_Widget* o = *a++;
			if (o == resizable()) continue;
			if (p[1]<q[1] && o->y()<=my+GRABAREA && o->y()+o->h()>=my-GRABAREA) {
				int t = mx - (o->x()+o->w());
				if (abs(t) < mindx) {
					sdx = t;
					mindx = abs(t);
					oldx = p[1];
				}
			}
			if (p[3]<q[3] && o->x()<=mx+GRABAREA && o->x()+o->w()>=mx-GRABAREA) {
				int t = my - (o->y()+o->h());
				if (abs(t) < mindy) {
					sdy = t;
					mindy = abs(t);
					oldy = p[3];
				}
			}
		}
		sdrag = 0; sx = sy = -1;
		if (mindx <= GRABAREA) {sdrag = DRAGH; sx = oldx;}
		if (mindy <= GRABAREA) {sdrag |= DRAGV; sy = oldy;}
		set_cursor(this, cursors[sdrag]);
		if (sdrag) return 1;
		return Fl_Group::handle(event);
	}

	case FL_LEAVE:
		set_cursor(this, FL_CURSOR_DEFAULT);
		break;

	case FL_DRAG:
		// This is necessary if CONSOLIDATE_MOTION in Fl_x.cxx is turned off:
		if (damage()) return 1; // don't fall behind
	case FL_RELEASE: {
		if (!sdrag) return 0; // should not happen
		Fl_Widget* r = resizable(); if (!r) r = this;
		int newx;
		if (sdrag&DRAGH) {
			newx = Fl::event_x()-sdx;
			if (newx < r->x()) newx = r->x();
			else if (newx >= r->x()+r->w()) newx = r->x()+r->w();
		} else
			newx = sx;
		int newy;
		if (sdrag&DRAGV) {
			newy = Fl::event_y()-sdy;
			if (newy < r->y()) newy = r->y();
			else if (newy >= r->y()+r->h()) newy = r->y()+r->h();
		} else
			newy = sy;
		position(sx,sy,newx,newy);
		if (event == FL_DRAG) set_changed();
		do_callback();
		return 1;}

	}

	return Fl_Group::handle(event);
}

