% Modify these two lines to reflect
% your account and password.
myaddress = 'opticaltable@gmail.com';
mypassword = 'ccdA126-Bilbo';

setpref('Internet','E_mail',myaddress);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',myaddress);
setpref('Internet','SMTP_Password',mypassword);

props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', ...
                  'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

sendmail({myaddress,'msalmanasif@gmail.com'}, 'Gmail Test', 'This is a test message.');
