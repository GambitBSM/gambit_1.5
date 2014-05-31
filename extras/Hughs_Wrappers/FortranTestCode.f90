program testcode
implicit none
double precision z, testfunction

call testsubroutine(1.0, 2.0)

z = testfunction(1.0, 2.0)

write (*,*) 'testfunction returned', z

end

subroutine testsubroutine (x, y)
real x, y

write (*,*) 'This is testsubroutine', x, y

end

double precision function testfunction(x, y)
real x, y

write (*,*) 'This is testfunction', x, y

testfunction = x + y

return

end

