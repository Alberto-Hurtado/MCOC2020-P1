import datetime as dt

utc_EOF_format= "%Y-%m-%dT%H:%M:%S.%f"
#t1=dt.datetime.strptime("2020-08-07T22:59:42.000000",utc_EOF_format)
t1=dt.datetime.strptime("2018-08-14T22:59:42.000000",utc_EOF_format) #profesor
#t2=dt.datetime.strptime("2020-08-09T00:59:42.000000",utc_EOF_format)
t2=dt.datetime.strptime("2018-08-16T00:59:42.000000",utc_EOF_format) #profesor

intervalo = t2-t1
intervalo_en_segundos=intervalo.total_seconds()