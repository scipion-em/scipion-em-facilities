# This is a template for the auxiliary file that contains
# the usernames, password and paths used to connecto to influx
# (influx section)
# and to transfer images between computers (paramiko section)
# The usernames, passwords, keyfilepath and keytype has been encrypted
# using the function enCrypt (see below)
# this encryption is weak but at least will stop casual users

# influx: information needed to acces to the "host"
# running influxdb. If you are not encrypting your
# communications set ssl = False
usernameInflux = 'aW5mbHV4dXNlcm5hbWU='
passwordInflux = 'aW5mbHV4cGFzc3dk'
dataBase = 'scipion'
hostinflux = 'influx-server.cnb.csic.es'
port = 8086
ssl = True
verify_ssl = False
timeZone = "Europe/Madrid"

# paramiko,  is a ssh client for python we use it to implement
# sftp and transfer images from scipion host to grafana host
# authentication is performed using username and a private key.
# The path to the private l¡key (keyfilepath) is encrypted and should be similar to
# '/home/transferusername/.ssh/id_rsa' and the keyfiletype (also encrypted)
# should be either "RSA" or "DSA"
# Remember to add the PUBLIC key to the authorized_host file in hostparamiko
usernameParamiko = 'dXNlcm5hbWVQYXJhbWlrbw=='
passwordParamiko = None,
keyfilepath = 'L2hvbWUvcm9iZXJ0by8uc3NoL2lkX3JzYQ=='
keyfiletype = 'UlNB'
remote_path = '/home/scipionbox/public_html/'
hostparamiko = 'paramiko-erver.cnb.csic.es'
