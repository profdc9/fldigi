// ----------------------------------------------------------------------------
//
// flxmlrpc Copyright (c) 2015 by W1HKJ, Dave Freese <iam_w1hkj@w1hkj.com>
//    
// XmlRpc++ Copyright (c) 2002-2008 by Chris Morley
//
// This file is part of fldigi
//
// flxmlrpc is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ----------------------------------------------------------------------------

#include <config.h>

#include "XmlRpcSource.h"
#include "XmlRpcUtil.h"

#if defined(_WINDOWS)
#undef _WINSOCKAPI_
# include <stdio.h>
# include <winsock2.h>
#else
extern "C" {
# include <unistd.h>
# include <stdio.h>
}
#endif

#if USE_OPENSSL_XMLRPC
# include <openssl/crypto.h>
# include <openssl/x509.h>
# include <openssl/pem.h>
# include <openssl/ssl.h>
# include <openssl/err.h>
#endif

namespace XmlRpc {

  struct SslProxy
  {
#if USE_OPENSSL_XMLRPC
    SSL_CTX* _ssl_ctx;
    SSL_METHOD* _ssl_meth;
    SSL* _ssl_ssl;
#endif
  };


  XmlRpcSource::XmlRpcSource(XmlRpcSocket::Socket fd /*= -1*/, bool deleteOnClose /*= false*/) 
    : _fd(fd)
    , _deleteOnClose(deleteOnClose)
    , _keepOpen(false)
    , _sslEnabled(false)
    , _ssl(0)
  {
  }

  XmlRpcSource::~XmlRpcSource()
  {
    delete _ssl;
  }


  void
  XmlRpcSource::setSslEnabled(bool b /*=true*/)
  {
#if USE_OPENSSL_XMLRPC
    _sslEnabled = b;
#endif
  }


  bool
  XmlRpcSource::doConnect()
  {
#if USE_OPENSSL_XMLRPC
    // Perform SSL if needed
    if (_sslEnabled)
    {
      _ssl = new SslProxy;

      SSLeay_add_ssl_algorithms();
      _ssl->_ssl_meth = const_cast<SSL_METHOD*>(SSLv23_client_method());
      SSL_load_error_strings();
      _ssl->_ssl_ctx = SSL_CTX_new(_ssl->_ssl_meth);
      _ssl->_ssl_ssl = SSL_new(_ssl->_ssl_ctx);
      SSL_set_fd(_ssl->_ssl_ssl, _fd);

      return SSL_connect(_ssl->_ssl_ssl) == 1;
    }
#endif
    return true;
  }

  // Read available text from the specified socket. Returns false on error.
  bool 
  XmlRpcSource::nbRead(std::string& s, bool *eof)
  {
    const int READ_SIZE = 4096;   // Number of bytes to attempt to read at a time
    char readBuf[READ_SIZE];

    bool wouldBlock = false;
    *eof = false;

    while ( ! wouldBlock && ! *eof)
    {
      int n;

#if USE_OPENSSL_XMLRPC
      // Perform SSL if needed
      if (_ssl && _ssl->_ssl_ssl)
      {
        n = SSL_read(_ssl->_ssl_ssl, readBuf, READ_SIZE-1);
      }
      else
#endif
#if defined(_WINDOWS)
      n = recv(_fd, readBuf, READ_SIZE-1, 0);
#else
      n = read(_fd, readBuf, READ_SIZE-1);
#endif
      XmlRpcUtil::log(5, "XmlRpcSocket::nbRead: read/recv returned %d.", n);

      if (n > 0) {
        readBuf[n] = 0;
        s.append(readBuf, n);
      } else if (n == 0) {
        *eof = true;
      } else if (XmlRpcSocket::nonFatalError()) {
        wouldBlock = true;
      } else {
        return false;   // Error
      }
    }
    return true;
  }


  // Write text to the socket. Returns false on error.
  bool 
  XmlRpcSource::nbWrite(std::string const& s, int *bytesSoFar)
  {
    int nToWrite = int(s.length()) - *bytesSoFar;
    const char *sp = s.c_str() + *bytesSoFar;
    bool wouldBlock = false;

    while ( nToWrite > 0 && ! wouldBlock )
    {
      int n;
#if USE_OPENSSL_XMLRPC
      // Perform SSL if needed
      if (_ssl && _ssl->_ssl_ssl)
      {
        n = SSL_write(_ssl->_ssl_ssl, sp, nToWrite);
      }
      else
#endif
#if defined(_WINDOWS)
      n = send(_fd, sp, nToWrite, 0);
#else
      n = write(_fd, sp, nToWrite);
#endif

      XmlRpcUtil::log(5, "XmlRpcSocket::nbWrite: send/write returned %d.", n);

      if (n > 0)
      {
        sp += n;
        *bytesSoFar += n;
        nToWrite -= n;
      }
      else if (XmlRpcSocket::nonFatalError())
      {
        wouldBlock = true;
      }
      else
      {
        return false;   // Error
      }
    }
    return true;
  }

  void
  XmlRpcSource::close()
  {
#if USE_OPENSSL_XMLRPC
    if (_ssl && _ssl->_ssl_ssl)
    {
      SSL_shutdown(_ssl->_ssl_ssl);

      // Should close be called here ? ...

      SSL_free(_ssl->_ssl_ssl);
      SSL_CTX_free(_ssl->_ssl_ctx);

      delete _ssl;
      _ssl = 0;
    }
#endif
    if ( (int)_fd != -1)
    {
      XmlRpcUtil::log(2,"XmlRpcSource::close: closing socket %d.", _fd);
      XmlRpcSocket::close(_fd);
      _fd = -1;
    }
    if (_deleteOnClose)
    {
      XmlRpcUtil::log(2,"XmlRpcSource::close: deleting this");
      _deleteOnClose = false;
      delete this;
    }
  }

} // namespace XmlRpc
