class PageController < ApplicationController
  before_filter :check_user
  
  
  def index
  end
  private
  
  def check_user
    if current_user.nil?
      redirect_to new_user_session_url
    end
  end
end
