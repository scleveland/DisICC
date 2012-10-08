# Sensors Demo
# Copyright (c) 2009 Montana State University
#
# FILE: voeis_controller.rb
# Controller for projects
class MonopredController < ApplicationController
  
  def index
    today = Time.now
    thirty_days = today + 30*24*60*60
    sixty_days = today + 60*24*60*60
    ninety_days = today + 90*24*60*60
    @thirty_day_items = Item.all(:review.gt => today, 
                                  :review.lt => thirty_days, 
                                  :order => [:review.desc])
    @sixty_day_items = Item.all(:review.gte => thirty_days, 
                                  :review.lt => sixty_days, 
                                  :order => [:review.desc])
    @ninety_day_items = Item.all(:review.gte => sixty_days, 
                                  :review.lt => ninety_days, 
                                  :order => [:review.desc])                                  
    @past_due = Item.all(:review.lte => today,:order => [:review.desc])
    
    respond_to do |format|
      format.html
    end
  end
end