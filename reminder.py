import smtplib
from email.mime.text import MIMEText
import sys,os

if len(sys.argv)< 2: 
  print("Usage: python reminder.py \"title\" \"content\"")
  exit(1)
 
#SMTP server provided by tencent_mail
mail_host = 'smtp.qq.com'
#Server port
port = 465
send_by = '1091325386@qq.com' # mail ID
password = 'hjcmhomuxxcojbcd' # STMP authorization code of QQ-mail
send_to = 'nyx233zmbxa@163.com' # send-to
def send_email(title,content,):
    #create MIMEText class, from class splain
      message = MIMEText(content,'plain','utf-8')
      message["From"] = send_by
      message['To'] = send_to
      message['Subject'] = title
      try:
          #mind 3rd parameter
          smpt = smtplib.SMTP_SSL(mail_host, port, 'utf-8')
          smpt.login(send_by,password)
          smpt.sendmail(send_by, send_to,message.as_string())
          print("Send email successfully!")
      except:
          print("Failed to send an email...")
def main():
    title = sys.argv[1]  # 1st input for title
    content = sys.argv[2]  # 2nd input for content
    send_email(title,content)
if __name__ == "__main__":
    main()