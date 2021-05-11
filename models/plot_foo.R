gg_fit_sampling <- function(data, name, lag_decrease = NULL, lag_stayhome = NULL) {
  p <- data %>% 
    ggplot() + 
    geom_point(aes(x=date, y=y)) + 
    geom_line(aes(x=date, y=fit_med), 
              col = "blue") + 
    geom_ribbon(aes(x=date, ymin=fit_lo, ymax=fit_hi), 
                alpha= 0.1, fill = "blue") + 
    labs(title = name, 
         x = "", y = "") +
    xlim(as.Date("2020-02-15"), as.Date("2020-06-30")) 

  if(!is.null(lag_decrease)) {
    p <- p + 
    geom_vline(aes(xintercept = decrease_50_total_visiting), color = "purple") + 
    geom_vline(aes(xintercept = decrease_50_total_visiting + lag_decrease), linetype="dotted", color = "purple")
  }

  if(!is.null(lag_stayhome)) {
    p <- p + 
    geom_vline(aes(xintercept = stayhome), color = "blue") + 
    geom_vline(aes(xintercept = stayhome + lag_stayhome), linetype="dotted", color = "blue")
  }

  p

}

gg_intrv_sampling <- function(data, name, lag_decrease = NULL, lag_stayhome = NULL) {
  p <- data %>% 
    ggplot() + 
    geom_point(aes(x=date, y=y)) + 
    geom_line(aes(x=date, y=fit_med), 
              col = "blue") + 
    geom_ribbon(aes(x=date, ymin=fit_lo, ymax=fit_hi), 
                alpha= 0.1, fill = "blue") + 
    geom_line(aes(x=date, y=ctr_med), 
              col = "red") + 
    geom_ribbon(aes(x=date, ymin=ctr_lo, ymax=ctr_hi), 
                alpha= 0.1, fill = "red") + 
    labs(title = name, 
         x = "", y = "") + 
    xlim(as.Date("2020-03-01"), as.Date("2020-04-30"))

  if(!is.null(lag_decrease)) {
    p <- p + 
    geom_vline(aes(xintercept = decrease_50_total_visiting), color = "purple") + 
    geom_vline(aes(xintercept = decrease_50_total_visiting + lag_decrease), linetype="dotted", color = "purple")
  }

  if(!is.null(lag_stayhome)) {
    p <- p + 
    geom_vline(aes(xintercept = stayhome), color = "blue") + 
    geom_vline(aes(xintercept = stayhome + lag_stayhome), linetype="dotted", color = "blue")
  }

  p

}

gg_intrv_effect <- function(data, name, lag_decrease = NULL, lag_stayhome = NULL) {
  p <- data %>% 
    ggplot() + 
    geom_point(aes(x=date, y=y_eff)) + 
    geom_line(aes(x=date, y=fit_med_eff), 
              col = "blue") + 
    geom_ribbon(aes(x=date, ymin=fit_lo_eff, ymax=fit_hi_eff), 
                alpha= 0.1, fill = "blue") + 
    geom_line(aes(x=date, y=ctr_med_eff), 
              col = "red") + 
    geom_ribbon(aes(x=date, ymin=ctr_lo_eff, ymax=ctr_hi_eff), 
                alpha= 0.1, fill = "red") + 
    labs(title = name, 
         x = "", y = "") + 
    xlim(as.Date("2020-03-01"), as.Date("2020-04-30"))

  if(!is.null(lag_decrease)) {
    p <- p + 
    geom_vline(aes(xintercept = decrease_50_total_visiting), color = "purple") + 
    geom_vline(aes(xintercept = decrease_50_total_visiting + lag_decrease), linetype="dotted", color = "purple")
  }

  if(!is.null(lag_stayhome)) {
    p <- p + 
    geom_vline(aes(xintercept = stayhome), color = "blue") + 
    geom_vline(aes(xintercept = stayhome + lag_stayhome), linetype="dotted", color = "blue")
  }

  p

}

gg_days_btwn_sampling <- function(data, name, up, down, lag_decrease = NULL, lag = NULL) {
  p <- data %>% 
    ggplot() + 
    geom_point(aes(x=date, y=y)) + 
    geom_line(aes(x=date, y=fit_med), 
              col = "blue") + 
    geom_ribbon(aes(x=date, ymin=fit_lo, ymax=fit_hi), 
                alpha= 0.1, fill = "blue") + 
      geom_line(aes(x=date, y=ctr1_med), 
              col = "red") + 
      geom_ribbon(aes(x=date, ymin=ctr1_lo, ymax=ctr1_hi), 
                alpha= 0.1, fill = "red") + 
       geom_line(aes(x=date, y=ctr3_med), 
              col = "green") + 
      geom_ribbon(aes(x=date, ymin=ctr3_lo, ymax=ctr3_hi), 
                alpha= 0.1, fill = "green") +  
      labs(title = name, x = "", y = "") + 
      xlim(as.Date("2020-03-01"), as.Date("2020-04-30"))
    
    if(!is.null(lag_decrease)) {
      p <- p + 
      geom_vline(aes(xintercept = decrease_50_total_visiting), color = "purple") + 
      geom_vline(aes(xintercept = decrease_50_total_visiting + lag_decrease), linetype="dotted", color = "purple") + 
      geom_vline(aes(xintercept = decrease_50_total_visiting + lag_decrease + down), linetype="dotted", color = "green") + 
      geom_vline(aes(xintercept = decrease_50_total_visiting + lag_decrease + up), linetype="dotted", color = "red")
    }
    
    if(!is.null(lag)) {
      p <- p + 
      geom_vline(aes(xintercept = stayhome), color = "blue") + 
      geom_vline(aes(xintercept = stayhome + lag), linetype="dotted", color = "blue") + 
      geom_vline(aes(xintercept = stayhome + lag + down), linetype="dotted", color = "green") + 
      geom_vline(aes(xintercept = stayhome + lag + up), linetype="dotted", color = "red")
    }

  p

}

gg_days_btwn_effect <- function(data, name, up, down, lag_decrease = NULL, lag = NULL) {
  p <- data %>% 
    ggplot() + 
    geom_point(aes(x=date, y=y_eff)) + 
    geom_line(aes(x=date, y=fit_med_eff), 
              col = "blue") + 
    geom_ribbon(aes(x=date, ymin=fit_lo_eff, ymax=fit_hi_eff), 
                alpha= 0.1, fill = "blue") + 
      geom_line(aes(x=date, y=ctr1_med_eff), 
              col = "red") + 
      geom_ribbon(aes(x=date, ymin=ctr1_lo_eff, ymax=ctr1_hi_eff), 
                alpha= 0.1, fill = "red") + 
       geom_line(aes(x=date, y=ctr3_med_eff), 
              col = "green") + 
      geom_ribbon(aes(x=date, ymin=ctr3_lo_eff, ymax=ctr3_hi_eff), 
                alpha= 0.1, fill = "green") +  
      labs(title = name, x = "", y = "") + 
      xlim(as.Date("2020-03-01"), as.Date("2020-04-30"))
    
    if(!is.null(lag_decrease)) {
      p <- p + 
      geom_vline(aes(xintercept = decrease_50_total_visiting), color = "purple") + 
      geom_vline(aes(xintercept = decrease_50_total_visiting + lag_decrease), linetype="dotted", color = "purple") + 
      geom_vline(aes(xintercept = decrease_50_total_visiting + lag_decrease + down), linetype="dotted", color = "green") + 
      geom_vline(aes(xintercept = decrease_50_total_visiting + lag_decrease + up), linetype="dotted", color = "red")
    }
    
    if(!is.null(lag)) {
      p <- p + 
      geom_vline(aes(xintercept = stayhome), color = "blue") + 
      geom_vline(aes(xintercept = stayhome + lag), linetype="dotted", color = "blue") + 
      geom_vline(aes(xintercept = stayhome + lag + down), linetype="dotted", color = "green") + 
      geom_vline(aes(xintercept = stayhome + lag + up), linetype="dotted", color = "red")
    }

  p

}

gg_covar_sampling <- function(data, name, covar_val, lag) {
  p <- data %>% 
    ggplot() + 
    geom_point(aes(x=date, y=y)) + 
    geom_line(aes(x=date, y=fit_med), 
              col = "blue") + 
    geom_ribbon(aes(x=date, ymin=fit_lo, ymax=fit_hi), 
                alpha= 0.1, fill = "blue") + 
    geom_line(aes(x=date, y=ctr1_med), 
              col = "red") + 
    geom_ribbon(aes(x=date, ymin=ctr1_lo, ymax=ctr1_hi), 
                alpha= 0.1, fill = "red") + 
    geom_line(aes(x=date, y=ctr3_med), 
              col = "green") + 
    geom_ribbon(aes(x=date, ymin=ctr3_lo, ymax=ctr3_hi), 
                alpha= 0.1, fill = "green") +  
    geom_vline(aes(xintercept = stayhome), color = "blue") + 
    geom_vline(aes(xintercept = stayhome + lag), linetype="dotted", color = "blue") +
    labs(title = paste(name, covar_val), 
         x = "", y = "") + 
    xlim(as.Date("2020-03-01"), as.Date("2020-04-30"))
  p
}

gg_compare_effect <- function(data, name) {
  p <- data %>% 
    ggplot() + 
    geom_point(aes(x=date, y=y_eff)) + 
    geom_line(aes(x=date, y=fit_med_eff_stayhome), 
              col = "blue") + 
    geom_line(aes(x=date, y=fit_med_eff_decrease), 
              col = "purple") +
    geom_line(aes(x=date, y=ctr_med_eff_stayhome), 
              col = "red") + 
    geom_line(aes(x=date, y=ctr_med_eff_decrease), 
              col = "orange") +
    labs(title = name, 
         x = "", y = "") + 
    xlim(as.Date("2020-03-01"), as.Date("2020-04-30")) + 
      geom_vline(aes(xintercept = decrease_50_total_visiting), color = "purple") + 
      geom_vline(aes(xintercept = decrease_50_total_visiting + 12), linetype="dotted", color = "purple") + 
      geom_vline(aes(xintercept = stayhome), color = "blue") + 
      geom_vline(aes(xintercept = stayhome + 12), linetype="dotted", color = "blue")
  p
}

gg_nchs_sampling <- function(nchs_data, n_, up = up, down = down) {
  
  p <- nchs_data %>% 
    ggplot() + 
    geom_point(aes(x=days_since_thresh, y=log_y)) + 
    geom_line(aes(x=days_since_thresh, y=fit_med), 
              col = "blue") + 
    geom_ribbon(aes(x=days_since_thresh, ymin=fit_lo, ymax=fit_hi), 
                alpha= 0.1, fill = "blue") + 
    geom_line(aes(x=days_since_thresh, y=ctr1_med), 
              col = "green") + 
    geom_ribbon(aes(x=days_since_thresh, ymin=ctr1_lo, ymax=ctr1_hi), 
                alpha= 0.1, fill = "green") + 
    geom_line(aes(x=days_since_thresh, y=ctr3_med), 
              col = "red") + 
    geom_ribbon(aes(x=days_since_thresh, ymin=ctr3_lo, ymax=ctr3_hi), 
                alpha= 0.1, fill = "red") +  
    #geom_vline(aes(xintercept = stayhome), color = "blue") + 
    #geom_vline(aes(xintercept = stayhome + 12), linetype="dotted", color = "blue") + 
    #geom_vline(aes(xintercept = stayhome + 12 + down), linetype="dotted", color = "green") + 
    #geom_vline(aes(xintercept = stayhome + 12 + up), linetype="dotted", color = "red") + 
    labs(title = n_, 
         x = "", y = "")
  p
}


gg_intrv_agg_sampling <- function(data, name, intrv_name, lag) {
  p <- data %>% 
    ggplot() + 
    geom_point(aes(x=days_since_thresh, y=y)) + 
    geom_line(aes(x=days_since_thresh, y=fit_med), 
              col = "blue") + 
    geom_ribbon(aes(x=days_since_thresh, ymin=fit_lo, ymax=fit_hi), 
                alpha= 0.1, fill = "blue") + 
    geom_line(aes(x=days_since_thresh, y=ctr_med), 
              col = "red") + 
    geom_ribbon(aes(x=days_since_thresh, ymin=ctr_lo, ymax=ctr_hi), 
                alpha= 0.1, fill = "red") + 
    xlim(-5, 40)
  
  if(intrv_name == "decrease") {
    p <- p + 
      geom_vline(aes(xintercept = days_btwn_decrease_thresh), color = "blue") + 
      geom_vline(aes(xintercept = days_btwn_decrease_thresh + lag), linetype="dotted", color = "blue") +
      labs(title = name, 
           x = "", y = "")
  }
  p
}
