/* MLPACK 0.2
 *
 * Copyright (c) 2008, 2009 Alexander Gray,
 *                          Garry Boyer,
 *                          Ryan Riegel,
 *                          Nikolaos Vasiloglou,
 *                          Dongryeol Lee,
 *                          Chip Mappus, 
 *                          Nishant Mehta,
 *                          Hua Ouyang,
 *                          Parikshit Ram,
 *                          Long Tran,
 *                          Wee Chin Wong
 *
 * Copyright (c) 2008, 2009 Georgia Institute of Technology
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 */
/**
 * @file common.c
 *
 * Implementations for bare-necessities FASTlib programming in C.
 */

#include "common.h"

#include <stdarg.h>
#include <unistd.h>

int segfault_on_abort = 0;
int abort_on_nonfatal = 0;
int pause_on_nonfatal = 0;
int print_notify_locs = 0;

char fl_msg_marker[] = {'X', '!', '*', '.'};
const char *fl_msg_color[] =
  {ANSI_HRED, ANSI_YELLOW, ANSI_HGREEN, ANSI_HBLUE};


int yet_another_isnan(double x) {
    return(x != x);
}

#ifdef WIN32
char FMT_LONG_LONG[] = "%I64d";
char FMT_ULONG_LONG[] = "%I64u";
#else
char FMT_LONG_LONG[] = "%lld";
char FMT_ULONG_LONG[] = "%llu";
#endif


void fl_abort(void)
{
  if (segfault_on_abort) {
    //fflush(NULL);
    //*(int*)NULL = 0;
  }
  //abort();
}

void fl_pause(void)
{
  if (isatty(0)) {
    char c;

    //fprintf(stderr, ANSI_HBLACK"Press Return to continue..."ANSI_CLEAR);
    //fflush(stderr);

    while ((c = getchar()) != EOF && c != '\n');
  }
}

void fl_print_msg_header(char marker, const char *color)
{
  //fprintf(stderr, "%s[%c]%s ", color, marker, ANSI_CLEAR);
}

void fl_print_msg_loc(const char *file, const char *func, int line)
{
  //const char *prev = file;
  const char *last = file;

  /* Finds the file and containing directory, if it exists. */
  while ((file = strchr(last, '/'))) {
    //prev = last;
    last = file + 1;
  }

  //fprintf(stderr, "%s:%s:%d: ", prev, func, line);
}

void fl_print_fatal_msg(const char *file, const char *func, int line,
                        const char *format, ...)
{
  va_list vl;

  fl_print_msg_header(fl_msg_marker[FL_MSG_FATAL],
		      fl_msg_color[FL_MSG_FATAL]);
  fl_print_msg_loc(file, func, line);

  va_start(vl, format);
  //vfprintf(stderr, format, vl);
  va_end(vl);

  //fprintf(stderr, "\n");

  fl_abort();
}

void fl_print_msg(const char *file, const char *func, int line,
                  fl_msg_t msg_type, const char *format, ...)
{
  va_list vl;

  fl_print_msg_header(fl_msg_marker[msg_type], fl_msg_color[msg_type]);
  if (msg_type < FL_MSG_NOTIFY_STAR || print_notify_locs) {
    fl_print_msg_loc(file, func, line);
  }

  va_start(vl, format);
  //vfprintf(stderr, format, vl);
  va_end(vl);

  //fprintf(stderr, "\n");

  if (msg_type < FL_MSG_NOTIFY_STAR) {
    if (msg_type == FL_MSG_FATAL || abort_on_nonfatal) {
      fl_abort();
    } else if (pause_on_nonfatal) {
      fl_pause();
    }
  }
}

void fl_print_progress(const char *desc, int prec)
{
  /*
  const int BAR_LEN = 50;
  static int prev_prec = -1;
  static const char *prev_desc = NULL;

  if (isatty(0)) {
    if (unlikely(prec != prev_prec || desc != prev_desc)) {
      char buf[BAR_LEN + 1];
      int pos = prec * BAR_LEN / 100;
      int i = 0;

      pos = unlikely(pos > BAR_LEN) ? BAR_LEN : unlikely(pos < 0) ? 0 : pos;

      for (; i < pos; ++i) {
        buf[i] = '#';
      }
      for (; i < BAR_LEN; ++i) {
        buf[i] = '.';
      }
      buf[i] = '\0';

      //fprintf(stderr, "\r"ANSI_BLUE"[%s] %d%% %s"ANSI_CLEAR"\r",
      //        buf, prec, desc);

      prev_prec = prec;
      prev_desc = desc;
    }
  }
   */
}

void hex_to_stream(FILE *stream, const char *src, const char *ok_char)
{
  char c;

  while ((c = *src++)) {
    if (isalnum(c) || strchr(ok_char, c)) {
      putc(c, stream);
    } else {
      fprintf(stream, "%%%02X", (unsigned)c);
    }
  }
}

char *hex_to_string(char *dest, const char *src, const char *ok_char)
{
  char c;

  while ((c = *src++)) {
    if (isalnum(c) || strchr(ok_char, c)) {
      *dest++ = c;
    } else {
      sprintf(dest, "%%%02X", (unsigned)c);
      dest += 3;
    }
  }

  *dest = '\0';
  return dest;
}

static int unhex(char c)
{
  if(c >= 'A' && c <= 'F') {
    return c-'A'+10;
  } else if(c >= 'a' && c <= 'f') {
    return c-'a'+10;
  } else {
    return c-'0';
  }
}

char *unhex_in_place(char *str)
{
  char *dest = strchr(str, '%');

  if (dest) {
    str = dest;
    while (*str) {
      if (*str == '%' && isxdigit(str[1]) && isxdigit(str[2])) {
        //printf("[%s:%d] cannot portably compute sscanf %%2hhx\n", __FILE__, __LINE__);
        //abort();
        *dest = unhex(str[1])*16 + unhex(str[2]);
        ++dest;
    	//sscanf(str + 1, FMT_HEX_CHAR, dest++);
	str += 3;
      } else {
	*dest++ = *str++;
      }
    }

    *dest = '\0';
    return dest;
  } else {
    return str + strlen(str);
  }
}
