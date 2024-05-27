remove_noise <- function(table){

  for (i in seq(1, length(table[, 1]), by=10)){       # Looping over every tenth entry in the table

    y = sum(table[i:(i+1000), 2]) / 1000              # Taking the mean of the next 1000 observations

    if(y > 39.99){                                    # If the mean is over 39.99 we break out of this loop
      break
    }
  }

  for (j in seq(length(table[, 1]), 1, by=-10)) {       # Looping over every tenth entry in the table, but this time we start at the end

    y = sum(table[(j-1000):j, 2]) / 1000                # Here we agian find the mean of the next 1000 observations

    if(y > 39.99){                                      # If the mean is over 39.99 we break out of this loop
      break
    }
  }

  return(subset(table, Time >= table[i, 1] & Time <= table[j, 1]))    # We return the index both loops reached before the mean was above 39.99
}



remove_repetition <- function(x, file_name){                               # takes case number and file name of excel file

  case = read_excel(file_name, sheet = "Målinger", guess_max = 100000)     # reads excel file with guess_max at 100000 to avoid boolean value if moniter has been running long before the data has been connected to the moniter.

  data_art <- VitalDBR::load_case(tname = 'SNUADC/ART', caseid = x)        # loads case from VitalDB to compare length of original data and data that has gone through the moniter

  data_art_wo_noice <- remove_noise(data_art)                              # removes noise from VitalDB data

  upsampled_data <- as.data.frame(seewave:::resamp(data_art_wo_noice$SNUADC.ART, 500, 1000, output='sample'))    # upsamples VitalDB data

  case <- case %>%
    select(MAP, HPI) %>%
    replace_with_na(replace = list(MAP = "mmHg"))                          # removes unit from first row

  case <- na.omit(case)                                                    # removes NAs

  n <- floor(nrow(upsampled_data)*0.001/20)                                # takes the numbers of rows of VitalDB data and times it by 0.001/20 to see the number of observations the data has when going from one sample pr. 0.002 sec to one sample pr. 20 sec

  if(n >= nrow(case)){                                                     # returns the entire moiter data if the moniter was stopped before the VitalDB data was done
    return(case)
  }

  else{
    return(case[1:n,])                                                    # returns the data with the length of the original data without repetition
  }

}




new_alarm_forward <- function(table, i) {

  # Code for when the current event is over in forward analysis, as defined when HPI falls below 85 again after an event
  # 'table' is the dataframe we want to loop through, and 'i' is the current index we have reached and observed an event

  for (k in seq(i, nrow(table))) {          # Looping over table from start value 'i' to end of the table
    if (table[k, 2] < 85) {                 # If the HPI in row 'k' is less than 85 we break from the loop
      break
    }
  }

  return(k)                                 # Return index k

}



new_alarm_backward <- function(table, i) {

  # Code for when the current event is over in backward analysis, as defined when ??? ???? ???? ??? again after an event
  # 'table' is the dataframe we want to loop through, and 'i' is the current index we have reached and observed an event

  for (k in seq(i, nrow(table))) {               # Looping over table from start value 'i' to end of the table
    #if (table[k, 2] < 85) {                   #new  - If the HPI in row 'k' is less than 85 we break from the loop
    if (table[k, 1] > 70) {                    #old -  If the MAP in row 'k' is greater than 70 we break from the loop
      break
    }
  }
  return(k)                                      # Return index k
}




biggest_dif <- function(table, i) {

  # Finds the biggest difference between MAPs in the table within a two minute interval
  # 'table' is the dataframe we want to loop through, and 'i' is the current index we have reached and want the see if there is a medical intervention

  biggest = 0

  for (k in seq(0, 5)) {

    if (table[i+k, 1] < 65 & i > 3) {
      #print(paste0("k = ", k))
      if(table[i+k-2, 1] < 65 & table[i+k-1, 1] < 65 & table[i+k, 1] < 65){
        return(FALSE)
      }
      if(table[i+k-1, 1] < 65 & table[i+k, 1] < 65 & table[i+k+1, 1] < 65){
        return(FALSE)
      }
      if(table[i+k, 1] < 65 & table[i+k+1, 1] < 65 & table[i+k+2 , 1] < 65){
        return(FALSE)
      }
    }

    for (m in seq(k+1, 6)) {

      if (table[i+m, 1] - table[i+k, 1] >= biggest ) {           #there is a difference of 8 or higher between to MAP points
        biggest = table[i+m, 1] - table[i+k, 1]

      }

    }

  }
  if (biggest >= 8){
    return(TRUE)
  }
  else {
    return(FALSE)
  }
}




fb_analysis <- function(case, table, minutes, hpi=TRUE){

  # This is a function that does a forward and backward analysis
  # Here we get case with is the number of the case that are being analised, 'table' with the columns MAP and HPI and minutes which is the time for the interval we are looking at in the analysis

  # Makes sure that the table columns are of the type numerical
  table$MAP <- as.numeric(table$MAP)
  table$HPI <- as.numeric(table$HPI)

  obs = minutes*3            # We have a sample for every 20 seconds, so to get the number of observations in the interval we multiply the minutes by 3

  n = nrow(table[, 1])       # Defining the number of observations in the table
  i = 1                      # The starting point of our first while loop

  #Initialising the types of alarms we can get
  LA = 0          # Late Alarm
  OA = 0          # Ongoing Alarm
  MI = 0          # Medical Intervention
  TP_f = 0        # True Positive in the forward analysis
  TP_b = 0        # True Positive in the backward analysis
  FP = 0          # False Positive
  FN = 0          # False Negative

  #Forward Analysis
  while (i <=  n-2) {      # In the loop we go from 1 to n-2 to not go out of bounds, as we have statements where we check two indexes ahead

    x = obs                # We define x as obs to allow for smaller intervals, to aviod getting out of bound

    if (i > n-obs) {       # If we are too close to the end of the table, x will be the amount of observations there are left
      x = n - i
    }


    if ((table[i, 2] >= 85 & table[i+1, 2] >= 85 & hpi) | (table[i, 1] <= 72 & table[i+1, 1] <= 72 &! hpi)) {   # If the observation we are at and the next is ≥ 85 we have an alarm
      print(i)
      for (j in seq(1, x)) {                         # At the alarm we will then go up to 'x' observations ahead to see if there is an event

        if (n-i < 4) {                             # Here we will be too close the the end of the table to check any kind of event so we break
          print(n-i)
          j = j + 1
          break
        }

        # Late alarm
        # If we have three consecutive observations where MAP < 65 and we are within one minute from when the alarm was rised we count it as a late alarm
        if (table[i+j, 1] < 65 & table[i+j+1, 1] < 65 & table[i+j+2, 1] < 65 & j < 4) {
          LA = LA + 1
          i = new_alarm_forward(table, i+j+1)      # The function new_alarm_forward finds the new index, which is when the current hypotensive event is over
          break
        }


        # True positives
        # If we have three consecutive observations where MAP < 65 and we are more than one minute from when the alarm was rised we count it as a true positive alarm
        if (table[i+j, 1] < 65 & table[i+j+1, 1] < 65 & table[i+j+2, 1] < 65 & j > 3) {
          TP_f = TP_f + 1
          i = new_alarm_forward(table, i+j+1)      # The function new_alarm_forward finds the new index, which is when the current hypotensive event is over
          break
        }


        # Medical intervention
        # If there is a difference of 5 or more between the MAP of the index we are at now and the next one or there is a difference of eight or higher within two minutes we define it as a medical intervention
        # The first part of both statements are to make sure we do not go out of bounds
        if (((n-i-j) >= 2 & (table[i+j+1, 1] - table[i+j, 1]) >= 5) | ((n-i-j) >= 6 & biggest_dif(table, i+j))) {
          MI = MI + 1
          i = new_alarm_forward(table, i+j+1)        # The function new_alarm_forward finds the new index, which is when the current event is over
          break
        }

        # False Positive
        # If we are at the last index of the loop, we have not had any of the above mentioned events, but an alarm was raised, so it must be a false positive alarm
        if (j == x) {
          FP = FP + 1
          i = new_alarm_forward(table, i+1)        # The function new_alarm_forward finds the new index, which is when the current event is over
          break
        }

      }

    }


    # When we have cheked for the above events we add 1 to 'i' to go to next index of the table
    else{
      i = i+1
    }

  }


  #Backward Analysis
  print("fw done")
  i = 6                   # We start at 6 because we do not count late alarms in the backwards analysis, so there is no reason to go through the first 5 observations

  while (i <= n-2) {      # We go to n-2 to not go out of bounds, as we check MAP from 'i' and two ahead

    if (table[i, 1] < 65 & table[i+1, 1] < 65 & table[i+2, 1] < 65) {   #if we have three consecutive MAP values under 65 we have a hypotensive event
      print(i)
      x = obs         # We define x as obs to allow for smaller intervals, to aviod getting out of bound

      if (i < obs) {  # If we are too close to the start of the table, x will be the amount of observations there are ahead of the index
        x = i - 3     # We subtract 3 to avoid going out of bounds, when checking if it is an ongoing alarm
      }


      for (j in rev(seq(2, x))) {    # Starting the loop at x and endning at 2, so we can count an alarm that starts just before one minute to the event

        if ((table[i-j, 2] >= 85 & table[i-j-1, 2] >= 85 & j >= 3 & hpi) | (table[i-j, 1] <= 72 & table[i-j-1, 1] <= 72 & j >= 3 &! hpi)) {        # there is an alarm that starts one minute or more before the event

          # Ongoning Alarm
          # If the alarm started before the interval we are checking in and all observations of HPI in the interval are greater than or equal to 85 we count it as an ongoing alarm
          if (table[i-j-2, 2] >= 85 & j == obs & all(table[i-x:i, 2] >= 85) ) {
            OA = OA + 1
            i = new_alarm_backward(table, i)         # The function new_alarm_backward finds the new index, which is when the current event is over
            break
          }

          # True positive
          # If there is an alarm in the interval and it isn't an ongoing alarm it is a true positive
          else{
            TP_b = TP_b + 1
            i = new_alarm_backward(table, i)         # The function new_alarm_backward finds the new index, which is when the current event is over
            break
          }

        }
        # False negative
        # If we have gone through the time interval and there hasn't been an alarm we count it as a false negative
        if (j == 2) {
          FN = FN + 1
          i = new_alarm_backward(table, i)           # The function new_alarm_backward finds the new index, which is when the current event is over
          break
        }

      }
    }

    # When we have cheked for the above events we add 1 to 'i' to go to next index of the table
    else{
      i = i+1
    }

  }

  # We return a list of case number and the number of alarms that has occurd
  return(list(case, LA, OA, MI, TP_f, TP_b, FP, FN))
}





confusion_matrix <- function(directory, minutes, HPI){


  columns= c("Case","LA", "OA", "MI", "TP_f", "TP_b", "FP", "FN")

  confusion = data.frame(matrix(nrow=1, ncol = length(columns)))

  colnames(confusion) = columns

  #getting a list of all files in the directory
  all_files <-list.files(directory)

  # used grep to find occurrences of strings that start with "Case" and end in ".xlsx"
  strings <- grep("^Case.*\\.xlsx$", all_files, value = TRUE)

  #looping over each file found in the directory that has monitor data
  for (string in strings){
    case_id <- as.numeric(strapplyc(string, "Case(\\d+)", simplify = TRUE))      #getting the case number from the file-name

    table <- remove_repetition(case_id, string)                                  #removing potential repetition in the monitor data

    confusion <- rbind(confusion, fb_analysis(case_id, table, minutes, hpi=HPI)) #getting the confusion values from the current file and adding it to the confusion. matrix
  }

  confusion <- confusion[-1, ]

  return(confusion)

}


calculate_ppv_sens <- function(table){
  ppv <- sum(table$TP_f)/(sum(table$TP_f) + sum(table$FP))
  sens <- sum(table$TP_b)/(sum(table$TP_b) + sum(table$FN))

  Positive_Predictive_Value= c(ppv)
  Sensitivity = c(sens)

  result <- data.frame(Positive_Predictive_Value, Sensitivity)

  return(result)
}











