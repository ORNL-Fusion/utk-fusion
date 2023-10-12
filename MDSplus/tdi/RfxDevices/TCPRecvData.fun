public fun TCPRecvData(in _sock)
{
	_dim = 0;
	_n_frame = 0;
	_x_pixel = 0;
	_y_pixel = 0;

	if( TcpClient->RecvData(val(_sock), ref(_dim), val(4)) != 4 )
	   return ( zero(0, 0B) );

write(*, "Buffer dimension ", _dim);


	_out = zero(_dim, 0B);
/*
	if( TcpClient->RecvData(val(_sock), ref(_out), val(_dim)) != _dim )
	   return ( zero(0, 0B) );
*/

	_dim = TcpClient->RecvData(val(_sock), ref(_out), val(_dim));

write(*, "Buffer dimension read", _dim);


	return ( _out );
}