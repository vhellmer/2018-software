from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import Select
from selenium.webdriver.support.ui import WebDriverWait as wait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.action_chains import ActionChains
from selenium.common.exceptions import TimeoutException
from time import sleep
import csv

#Gets list of all iGEM teams from csv file
teamlist = []
with open('2018__team_list__2018-09-10.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
    	print(row)
    	if row[7] == 'Accepted':
    		teamlist.append(row[1])


#Opens chrome (MUST have selenium plugin in PATH or will cause an error)
driver = webdriver.Chrome()

#File for strains to be written to
txtfile = open("strainfile.txt", "w")

#Gets strain information from every team. Reloads page and handles a pop-up which appears for going into viewing mode.
#The sleep() commands are so information can load and so their sever doesn't get overloaded and fail (there are better ways to do this, but this works for now)
for team in teamlist:
		sleep(2)
		driver.get("http://2018.igem.org/Safety/Final_Safety_Form")
		select = Select(driver.find_element_by_id('team_list_dropdown'))
		s = driver.find_element_by_id('team_list_dropdown')
		s.click()
		sleep(1)

		#Trys to find team in dropdown menu and strain info, otherwise prints teamname if an exception is raised
		try:
			select.select_by_visible_text(str(team))
			sleep(2)
			alert = driver.switch_to.alert
			alert.accept()
			sleep(2)
			textarea = driver.find_element_by_xpath('//*[@id="formbody"]/fieldset[3]/div/div[1]/textarea').get_attribute("value")
			if(textarea != ""):
				txtfile.write(textarea + "\n\n")
		except:
			print(team + " not found in dropdown menu")

txtfile.close()
driver.close()

